#!/usr/bin/env python
###############################################################################

import sys

import fastcluster as fc
import numpy as np

import utils as utils
import blastutils as Butils
import dendrogramutils as Dutils

###############################################################################

hitScoreHitsFile = sys.argv[1]
exonInfoFile     = sys.argv[2]
dendrogramFile   = sys.argv[3]

###############################################################################
  # Read the input files
H = Butils.readBlastFile(hitScoreHitsFile, Butils.genomicblastfields)
exonInfo = utils.readExonInfo(exonInfoFile)

###############################################################################

  # Generate a distance matrix for each pair of hits

def distance(hit_i, hit_j):

  def dist(r1, r2):
    (s1, e1) = r1
    (s2, e2) = r2
    # <= do not allow same start as end end. (100, 200), (200, 300) NOT ALLOWED
    # <  allow same start as end. (100, 200), (200, 300) ALLOWED
    return max(0, (max(s1, s2) - min(e1, e2)))
  #edef


  # Refer to eq. 2 in the paper

  # CHECK FOR STRANDYNESS
  rbi = (hit_i.g1start, hit_i.g1end)
  rbj = (hit_j.g1start, hit_j.g1end)

  rgi = (hit_i.g2start, hit_i.g2end)
  rgj = (hit_j.g2start, hit_j.g2end)

  return dist(rbi, rbj) + dist(rgi, rgj)

#edef

D = []

for i in range(len(H)):
  for j in range(i+1, len(H)-1):
    D.append(distance(H[i], H[j]))
  #efor
#efor

###############################################################################

dendrogram = []
if len(D) > 0:
  dendrogram = fc.linkage(D, method='single')
#fi

tree = Dutils.prepareDendrogram(dendrogram)

Dutils.writeDendrogram(tree, dendrogramFile)

