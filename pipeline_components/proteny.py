#!/usr/bin/env python

import sys

import utils as utils
import blastutils as Butils
import dendrogramutils as Dutils
import proteny_core as pcore

###############################################################################
  # Collect the command line parameters

dendrogramHitsListFile = sys.argv[1]
exonInfoFile           = sys.argv[2]
pvalue                 = float(sys.argv[3])
cons_thresh            = float(sys.argv[4])
threads                = int(sys.argv[5])
outputFile             = sys.argv[6]

###############################################################################
  # Read the data

PS = pcore.ProtenyStructure(dendrogramHitsListFile, exonInfoFile, threads)

###############################################################################
  # Detect significant clusters

sigClusters = PS.getSigClusters(pvalue, cons_thresh)

###############################################################################

print(outputFile)
pcore.writeClusterInfo(PS, sigClusters, outputFile)
