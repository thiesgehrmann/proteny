###############################################################################

import errno
import os
import csv
import gzip
from collections import namedtuple

###############################################################################

blastfields = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"
blastfields_type = "str str float int int int int int int int float float int int"
augmentedblastfields = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen K"
augmentedblastfields_type = "str str float int int int int int int int float float int int float"
genomicblastfields = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen K g1chr g1start g1end g2chr g2start g2end"
genomicblastfields_type = "str str float int int int int int int int float float int int float str int int str int int"


###############################################################################

# https://stackoverflow.com/a/600612
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

###############################################################################

class BlastHitType(object):
  BlastHit = None
  def __init__(self, blastFields = blastfields):
    self.BlastHit = namedtuple("BlastHit", blastFields)
  def setFields(self, blastFields):
    self.BlastHit = namedtuple("BlastHit", blastFields)
  #edef
  def getFields(self):
    return self.BlastHit._fields

  def blastFieldTypeConversion(self, fieldName, value):
    typeMap = { "qseqid" : lambda x: str(x),
                "sseqid" : lambda x: str(x),
                "pident" : lambda x: float(x),
                "length" : lambda x: int(x),
                "mismatch" : lambda x: int(x),
                "gapopen": lambda x: int(x),
                "qstart": lambda x: int(x),
                "qend": lambda x: int(x),
                "sstart": lambda x: int(x),
                "send": lambda x: int(x),
                "evalue": lambda x: float(x),
                "bitscore": lambda x: float(x),
                "slen": lambda x: int(x),
                "qlen": lambda x: int(x),
                "K": lambda x: float(x),  # AUGMENTED
                "g1chr": lambda x: str(x), # GENOMIC
                "g1start" : lambda x: int(x), #GENOMIC
                "g1end": lambda x: int(x), #GENOMIC
                "g2chr": lambda x: str(x), #GENOMIC
                "g2start" : lambda x: int(x), #GENOMIC
                "g2end": lambda x: int(x) } #GENOMIC
    return typeMap[fieldName](value)
  #edef
#eclass

blastHitType = BlastHitType()

def parseBlastHit(row, blastType=blastHitType):
  typedRow = [ blastType.blastFieldTypeConversion(field, value) for (field, value) in zip(blastType.getFields(), row) ]
  return blastType.BlastHit(*typedRow)
#edef

###############################################################################

def blastHit2Row(h):
  return [str(x) for x in [h[i] for i in range(len(h._fields))]]
#edef

###############################################################################

def readBlastFile(filename, fields=blastfields):
  blastHitType.setFields(fields)
  hits    = []
  # Read the BLAST hits
  with open(filename, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t", quotechar="\"");
    for row in reader:
      hits.append(parseBlastHit(row))
    #efor
  #ewith
  return hits
#edef

###############################################################################

def writeBlastFile(H, filename):
  with open(filename, "w") as ofd_hits:
    for h in H:
      ofd_hits.write('\t'.join(blastHit2Row(h)) + '\n')
    #efor
  #ewith
#edef

###############################################################################

def blast_hits_overlap(hit_i, hit_j, logic=lambda x,y : x or y):
  # note about logic:
  # OR: overlapping in query OR subject
  # and: overlapping in query AND subject

  if ((hit_i.qseqid == hit_j.qseqid) and (hit_i.sseqid == hit_j.sseqid)):
    queryOverlap   = (max(hit_i.qstart, hit_j.qstart) <= min(hit_i.qend, hit_j.qend))
    subjectOverlap = (max(hit_i.sstart, hit_j.sstart) <= min(hit_i.send, hit_j.send))
    return logic(queryOverlap, subjectOverlap)
  elif ((hit_i.qseqid == hit_j.sseqid) or (hit_i.sseqid == hit_j.qseqid)):
    queryOverlap   = (max(hit_i.qstart, hit_j.qstart) <= min(hit_i.send, hit_j.send))
    subjectOverlap = (max(hit_i.sstart, hit_j.sstart) <= min(hit_i.qend, hit_j.qend))
    return logic(queryOverlap, subjectOverlap)
  else:
    return False
  #fi
#edef

def regions_overlap(r1, r2):
  (s1, e1) = r1
  (s2, e2) = r2
  # <= do not allow same start as end end. (100, 200), (200, 300) NOT ALLOWED
  # <  allow same start as end. (100, 200), (200, 300) ALLOWED
  return max(0, (min(e1, e2) - max(s1, s2)))
#edef

###############################################################################

def blast_hits_overlap_group_helper(hit_i, hit_j, allowed_overlap=5):
  same_query   = hit_i.qseqid == hit_j.qseqid
  same_subject = hit_i.sseqid == hit_j.sseqid

  query_overlap   = regions_overlap((hit_i.qstart, hit_i.qend), (hit_j.qstart, hit_j.qend)) > allowed_overlap
  subject_overlap = regions_overlap((hit_i.sstart, hit_i.send), (hit_j.sstart, hit_j.send)) > allowed_overlap

  return (same_query and query_overlap) or (same_subject and subject_overlap)
#wedef

def blast_hits_overlap_group(hits):
  return any([ blast_hits_overlap_group_helper(hits[i],hits[j]) for i in range(len(hits)-1) for j in range(i+1,len(hits)) ])
#edef

###############################################################################

def bestHitPerSequence(hits, func=lambda hit: hit.bitscore):
  bestHit = {}
  for hit in hits:
    if hit.qseqid not in bestHit:
      bestHit[hit.qseqid] = func(hit)
    else:
      if func(hit) > bestHit[hit.qseqid]:
        bestHit[hit.qseqid] = func(hit)
      #fi
    #fi

    if hit.sseqid not in bestHit:
      bestHit[hit.sseqid] = func(hit)
    else:
      if func(hit) > bestHit[hit.sseqid]:
        bestHit[hit.sseqid] = func(hit)
      #fi
    #fi
  #efor
  return bestHit
#edef

