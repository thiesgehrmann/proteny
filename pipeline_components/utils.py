###############################################################################

import errno
import os
import csv
import gzip
from collections import namedtuple

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

def loadFasta(fastaFile):

  F = {'': 0}

  current_seq = ""
  buffer_seq  = ""
  
  with (gzip.open(fastaFile, "r") if fastaFile[-2:] == "gz" else open(fastaFile, "r")) as fd:
    for line in fd:
      line = line.strip()
      if len(line) == 0:
        continue
      #fi
      if line[0] == '>':
        F[current_seq] = buffer_seq
        current_seq = line[1:].split(' ')[0]
        buffer_seq = ""
      else:
        buffer_seq = buffer_seq + line.strip()
      #fi
  #ewith
  F[current_seq] = buffer_seq
  F.pop("", None)
  return F
#edef

###############################################################################

def writeFasta(fasta, outFile, linelength=80):
  with open(outFile, "w") as ofd:
    for  (name, sequence) in fasta:
      ofd.write(">%s\n" % name)
      ofd.write("%s\n" % '\n'.join([sequence[i:i+linelength] for i in range(0, len(sequence), linelength)]))
    #efor
  #ewith
#edef

###############################################################################

import numpy as np
def region_covered(regions, size=None):
  if size is None:
    size = max([ r[1] for r in regions ])
  #fi
  R = np.zeros(size)
  for r in regions:
    R[r[0]-1:r[1]-1] = 1
  #efor
  return R
#efor

###############################################################################

def region_coverage(regions, size = None):
  if size is None:
    size = max([ r[1] for r in regions ])
  #fi
  R = np.zeros(size)
  for r in regions:
    R[r[0]-1:r[1]-1] += 1
  #efor
  return R
#edef

###############################################################################

def indexListBy(L, key=lambda x: x[0]):
  G = {}
  for item in L:
    k = key(item)
    if k not in G:
      G[k] = []
    #fi
    G[k].append(item)
  #efor
  return G
#edef


###############################################################################


def readMapping(mapping):
  import csv
  M = []
  with open(mapping, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t")
    for (g1, g2) in reader:
      M.append((g1,g2))
    #efor
  #ewith
  return M
#edef

###############################################################################

def readColumnFile(filename, columnNames, delimiter='\t', types=""):
  import csv
  L = []
  typeFunctions = { "str" : lambda x: str(x),
                    "int" : lambda x: int(x),
                    "float" : lambda x: float(x) }

  if types != "":
    types = [ typeFunctions[c] for c in types.split(" ") ]
  #fi

  lineType = namedtuple("lineType", columnNames)
  with open(filename, "r") as ifd:
    reader = csv.reader(ifd, delimiter=delimiter)
    for row in reader:
      if row[0][0] == '#':
        continue
      #fi
      if len(types) == len(row):
        row = [ tf(v) for (tf, v) in zip(types, row) ]
      #fi
      L.append(lineType(*row))
    #efor
  #ewith
  return L
#edef

###############################################################################

def readExonInfo(filename):
  return dict([ (e.exonid, e) for e in readColumnFile(filename, "exonid proteinid genome seqid start end strand", types="str str str str int int str") ])
#edef

def indexExonInfo(exonInfo):
  from intervaltree import Interval, IntervalTree

  genomeChromosomes = dict([])

  for exonID in exonInfo.keys():
    exon = exonInfo[exonID]

    if exon.genome not in genomeChromosomes:
      genomeChromosomes[exon.genome] = dict([])
    #fi
    if exon.seqid not in genomeChromosomes[exon.genome]:
      genomeChromosomes[exon.genome][exon.seqid] = IntervalTree()
    #fi

    if exon.start != exon.end:
      genomeChromosomes[exon.genome][exon.seqid][exon.start:exon.end] = exon.exonid
    #fi
  #efor

  return genomeChromosomes

#edef

###############################################################################
