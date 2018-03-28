import errno
import os
import csv
import gzip
from collections import namedtuple


import gzip

###############################################################################

gff3Entry = namedtuple("gff3Entry", "seqid, source, type, start, end, score, strand, phase, attr")

def parseGFF3entry(fields):

   (seqid, source, feature, start, end, score, strand, phase, attr) = fields
   attr = [ x.strip().split('=') for x in attr.split(";") ]
   attr = dict([ (a[0], a[1]) if len(a) > 1 else (a[0], None) for a in attr ])
   return gff3Entry(seqid, source, feature, int(start), int(end), score, strand, phase, attr)
#edef

###############################################################################

class GFF3(object):
  def __init__(self, entries):
    self.entries = entries
    self.seqids  = set([ e.seqid for e in self.entries])
    self.genes   = { e.attr["ID"] : e for e in self.entries if e.type.lower() == 'mrna' }
    self.interval = self.indexByInterval()
    #self.geneIndex = self.indexByTopLevel()
  #edef

  def indexByTopLevelID(self):
    topLevel = { e.attr["ID"]: [e] for e in self.entries if "Parent" not in e.attr }
    lowerLevel = { e.attr["ID"] : e.attr["Parent"] for e in self.entries if "Parent" in e.attr }
    for e in [ e for e in self.entries if "Parent" in e.attr ]:
      parent = e.attr["Parent"]
      while parent not in topLevel:
        if parent in lowerLevel:
          parent = lowerLevel[parent]
        else:
          print("%s NOT in lowerLevel!" % parent)
          break
        #fi
      #ewhile
      topLevel[parent].append(e)
    #efor
    return topLevel
  #edef


    # I use this code inside a conda environment, but it is also loaded in the regular snakemake file.
    # To prevent this from being a problem, I allow it to gracefully fail when the intervaltree library doesn't exist
  def indexByInterval(self):
    try: 
      from intervaltree import Interval, IntervalTree
    
      t = { seqid: IntervalTree() for seqid in self.seqids }
      for e in [ e for e in self.entries if e.type.lower() == "mrna"]:
        t[e.seqid][e.start:e.end] = e
      #efor
      return t
    except ImportError:
      return {}
  #edef

  def areTandem(self, id1, id2):
    gene1 = self.genes[id1]
    gene2 = self.genes[id2]
    if gene1.seqid != gene2.seqid:
      return False
    #fi

    inregion = set([ e[-1].attr["ID"] for e in self.interval[gene1.seqid][min(gene1.end,gene2.end):max(gene1.start,gene2.start)] ])
    if len(inregion - set([gene1.attr["ID"], gene2.attr["ID"]])) == 0:
      return True
    else:
      return False
    #fi
  #edef

  def areSameStrand(self, id1, id2):
    gene1 = self.genes[id1]
    gene2 = self.genes[id2]
    return gene1.strand == gene2.strand
  #edef


###############################################################################

def readGFF3File(filename):
  G = []
  with (gzip.open(filename, "rt") if filename[-2:] == "gz" else open(filename, "r")) as gffFile:
    gffReader = csv.reader(gffFile, delimiter="\t", quotechar='"')
    for row in gffReader:
      if len(row) != 9:
        continue
      #fi
      G.append(parseGFF3entry(row))
    #efor
  #ewith
  return GFF3(G)
#edef

###############################################################################
