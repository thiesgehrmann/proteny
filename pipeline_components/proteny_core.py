import math
import bisect
import multiprocessing

from scipy import stats;
import numpy as np

import utils as utils
import blastutils as Butils
import dendrogramutils as Dutils

###############################################################################

class ProtenyStructure(object):
  #exonInfo = None
  #D     = {}
  #bestHits = {}
  #testCounter = 0
  #alpha = 0.05
  #conservationThreshold = 2.0
  #self.permutations = {}

  #A      = None
  #meanA  = None
  #varA   = None
  #U1     = None
  #meanU1 = None
  #varU1  = None
  #U2     = None
  #meanU2 = None
  #varU2  = None

  #parPool = None

  ###############################################################################

  def __init__(self, dendrogramHitsListFile, exonInfoFile, threads):
    self.exonInfo = utils.readExonInfo(exonInfoFile)
    self.exonInfoTree = utils.indexExonInfo(self.exonInfo)

    dendrogramHitsList = utils.readColumnFile(dendrogramHitsListFile, "chrA chrB hitsfile dendrogramfile")
    self.A = []
    self.D = {}
    for dh in dendrogramHitsList:
      hits          = Butils.readBlastFile(dh.hitsfile, Butils.genomicblastfields)
      dendrogram    = Dutils.readDendrogram(dh.dendrogramfile)
      nodeTestInfos = [ None for n in dendrogram ]
      self.D[(dh.chrA, dh.chrB)] = { "hits" : hits,
                             "dendrogram" : dendrogram,
                             "nodetestinfo" : nodeTestInfos }
      self.A.extend([ h.K for h in hits])

    #efor

    self.bestHits = Butils.bestHitPerSequence([ hit for hits in [ self.D[chrpair]["hits"] for chrpair in self.D.keys() ] for hit in hits ], func=lambda x: x.K)

    self.U1 = []
    self.U2 = []
    for exon in self.exonInfo.values():
      if (exon.exonid not in self.bestHits):
        continue
      #fi

      if exon.genome == "genome_1":
        self.U1.append(self.bestHits[exon.exonid])
      else:
        self.U2.append(self.bestHits[exon.exonid])
      #fi
    #efor

    self.meanA  = np.mean(self.A)
    self.varA   = np.var(self.A)
    self.meanU1 = np.mean(self.U1)
    self.varU1  = np.var(self.U1)
    self.meanU2 = np.mean(self.U2)
    self.varU2  = np.var(self.U2)

    self.permutations = {}
    self.parPool = multiprocessing.Pool(threads)
    self.testCounter = 1

  #edef

  ###############################################################################

  def getSigClusters(self, alpha, conservationThreshold):

    self.alpha = alpha
    self.conservationThreshold = conservationThreshold

    sigClusts = set([ (chrpair, -1) for chrpair in self.D.keys() ])
    prevSigClusts = []
    
    while sigClusts != prevSigClusts:
      print("#############################\n\n\n\n\nDOING IT AGAIN\n\n\n\n\n################################")
      prevSigClusts = sigClusts
      newSigClusts = []
    
      for (chrpair, nodeID) in sigClusts:
        newSigClusts.extend(self.getSigClustersAtNode(chrpair, nodeID))
      #efor
      sigClusts = set(newSigClusts)
    #ewhile

    return sigClusts
  #edef

  ###############################################################################

  def getSigClustersAtNode(self, chrpair, startNodeID):
    nodeList = [ startNodeID ]
    sigClusters = []

    while len(nodeList) > 0:
      nodeID = nodeList.pop()

      node = self.getNode(chrpair, nodeID)

      self.testNode(chrpair, nodeID)
      self.testNode(chrpair, node.child1)
      self.testNode(chrpair, node.child2)

      nodeInfo = self.getNodeInfo(chrpair, nodeID)
      print("(%s, %s)[%d]-> S:%f C:%f T:%s Z:%f P:%f" % (chrpair[0], chrpair[1], nodeID, nodeInfo.score, nodeInfo.conservationScore, "CLT" if nodeInfo.cltTest else "PERM", nodeInfo.zscore, nodeInfo.pvalue))

        # If the node is significant
      if nodeInfo.isSignificant(self.alpha, self.testCounter):
        print("Node is significant")
        # We add a node ONLY if:
        #  * it is significant
        #  * children are less significant
        #  * children are eligible nodes (more than 1 gene, conservation threshold is satisfactory)
        child1NodeInfo = self.getNodeInfo(chrpair, node.child1)
        child2NodeInfo = self.getNodeInfo(chrpair, node.child2)

        zscore = nodeInfo.zscore
        if child1NodeInfo.zscore >= zscore:
          sigClusters.append((chrpair, node.child1))
        #fi

        if child2NodeInfo.zscore > zscore:
          sigClusters.append((chrpair, node.child2))
        #fi

        if zscore >= child1NodeInfo.zscore and zscore > child2NodeInfo.zscore:
          sigClusters.append((chrpair, nodeID))
        #fi
        
        # If it is not significant, then return the significant nodes of the children nodes
      else:
        if node.child1 != -1:
          nodeList.append(node.child1)
        #fi
        if node.child2 != -1:
          nodeList.append(node.child2)
      #fi
    #ewhile
    return sigClusters
      
  #edef

  ###############################################################################

  def testNode(self, chrpair, nodeID):
    if self.getNodeInfo(chrpair, nodeID) is None:
      (genes1, genes2) = self.getNodeGenes(chrpair, nodeID)
      nGenes1 = len(genes1)
      nGenes2 = len(genes2)

      zscore = 0.0
      pvalue = 1.0

        # If the node doesn't have at least two genes in one genome, ignore it!
      if (nGenes1 < 2) and (nGenes2 < 2):
        self.D[chrpair]["nodetestinfo"][nodeID] = NodeTestInfo(nGenes1, nGenes2, 0.0, 0, 0, 0, 0, True, zscore, pvalue, self.testCounter)
        return
      #fi

      nA, nU1, nU2, conservationScore, nodeScore = self.nodeScore(chrpair, nodeID)

        # If the conservationScore is too low
      if conservationScore < self.conservationThreshold:
        self.D[chrpair]["nodetestinfo"][nodeID] = NodeTestInfo(nGenes1, nGenes2, conservationScore, 0, 0, 0, 0, True, zscore, pvalue, self.testCounter)
        return
      #fi

      self.increaseTestCounter()

      cltTest = (nA >= 10) and (nU1 >= 10) and (nU2 >= 10)

      zscore = 0.0
      pvalue = 1.0
      if cltTest:
        zscore, pvalue = self.cltTest(nA, nU1, nU2, nodeScore)
      else:
        zscore, pvalue = self.permutationTest(nA, nU1, nU2, nodeScore)
      #fi

      self.D[chrpair]["nodetestinfo"][nodeID] = NodeTestInfo(nGenes1, nGenes2, conservationScore, nA, nU1, nU2, nodeScore, cltTest, zscore, pvalue, self.testCounter)
    else:
      nodeInfo = self.getNodeInfo(chrpair, nodeID)
        # If it is CLT, we don't need to do anything anymore.
        # But, if it is not CLT, we may need to do some more permutations.
        # ONLY IF IT WAS SIGNIFICANT the last time we tested!
      if not(nodeInfo.cltTest) and nodeInfo.isSignificant(nodeInfo.pvalue, nodeInfo.lastNTests):
        zscore, pvalue = self.permutationTest(nodeInfo.nA, nodeInfo.nU1, nodeInfo.nU2, nodeInfo.score)
        nodeInfo.setpvalue(pvalue)
        nodeInfo.setzscore(zscore)
        nodeInfo.setlastNTests(self.testCounter)
      #fi
    #fi
        
  #edef

  ###############################################################################

  def getNodeGenes(self, chrpair, nodeID):
    genes1 = []
    genes2 = []
    hits    = self.getNodeHits(chrpair, nodeID)
    for h in hits:
      genes1.append(self.exonInfo[h.qseqid].proteinid)
      genes2.append(self.exonInfo[h.sseqid].proteinid)
    #efor

    return (set(genes1), set(genes2))
  #edef

  ###############################################################################

  def getNodeRegions(self, chrpair, nodeID):
    nodeHits = self.getNodeHits(chrpair, nodeID)

    genome1Start = min([ h.g1start for h in nodeHits ])
    genome1End   = max([ h.g1end for h in nodeHits ])
    genome2Start = min([ h.g2start for h in nodeHits ])
    genome2End   = max([ h.g2end for h in nodeHits ])

    return ((genome1Start, genome1End), (genome2Start, genome2End))
  #edef

  ###############################################################################

  def getNodeHits(self, chrpair, nodeID):
    return [ self.D[chrpair]["hits"][hitID] for hitID in self.getNode(chrpair, nodeID).leaves ]
  #edef

  ###############################################################################

  def getNodeInfo(self, chrpair, nodeID):
    return self.D[chrpair]["nodetestinfo"][nodeID]
  #edef

  ###############################################################################

  def getNode(self, chrpair, nodeID):
    return self.D[chrpair]["dendrogram"][nodeID]
  #edef

  ###############################################################################

  def nodeScore(self, chrpair, nodeID):
    node = self.getNode(chrpair, nodeID)
    nodeInfo = self.getNodeInfo(chrpair, nodeID)
    if nodeInfo is not None:
      return (nodeInfo.nA, nodeInfo.nU1, nodeInfo.nU2, node.conservationScore, nodeInfo.score)
    #fi

    nodeHits = self.getNodeHits(chrpair, nodeID)
    ((genome1Start, genome1End), (genome2Start, genome2End)) = self.getNodeRegions(chrpair, nodeID)

    exonsInGenome1Region = [ r[-1] for r in self.exonInfoTree["genome_1"][chrpair[0]][min(genome1Start, genome1End):max(genome1Start, genome1End)] ]
    exonsInGenome2Region = [ r[-1] for r in self.exonInfoTree["genome_2"][chrpair[1]][min(genome2Start, genome2End):max(genome2Start, genome2End)] ]

    accountedExons          = set([ h.qseqid for h in nodeHits] + [ h.sseqid for h in nodeHits])
    unaccountedExonsGenome1 = [ e for e in exonsInGenome1Region if e not in accountedExons ]
    unaccountedExonsGenome2 = [ e for e in exonsInGenome2Region if e not in accountedExons ]
    unaccountedExonsGenome1Elsewhere = [ e for e in unaccountedExonsGenome1 if e in self.bestHits ]
    unaccountedExonsGenome2Elsewhere = [ e for e in unaccountedExonsGenome2 if e in self.bestHits ]

    sumK = 2 * sum([ hit.K for hit in nodeHits ])

    unaccountedK = sum([ self.bestHits[ue] for ue in unaccountedExonsGenome1 + unaccountedExonsGenome2 if ue in self.bestHits ])

    #print("%s,%s (%d) -> %f" % (chrpair[0], chrpair[1], nodeID, sumK - unaccountedK))
    score = sumK - unaccountedK

    nA = len(accountedExons)
    nU1 = len(unaccountedExonsGenome1)
    nU2 = len(unaccountedExonsGenome2)
    nU1EO = len(unaccountedExonsGenome1Elsewhere)
    nU2EO = len(unaccountedExonsGenome2Elsewhere)

    conservationScore = nA / float(nU1EO + nU2EO + 1)

    return (nA, nU1, nU2, conservationScore, score)
  #edef

  ###############################################################################

  def cltTest(self, nA, nU1, nU2, score):
    mu = 2*(self.meanA * nA) - (self.meanU1 * nU1) - (self.meanU2 * nU2);
    s2 = (4 * self.varA * nA) + (self.varU1 * nU1) + (self.varU2 * nU2);

    z = ( score - mu ) / math.sqrt(s2);
    p = 1 - stats.norm.cdf(z);

    #print("pvalue -> %f -> (%f, %f) -> %f"% (score, mu, s2, p))
 
    return (z,p)
  #edef

  ###############################################################################

  def permutationTest(self, nA, nU1, nU2, score):
    import math
    if (nA, nU1, nU2) not in self.permutations:
      self.permutations[(nA,nU1,nU2)] = [ 0, [] ] # [ total, highest ]
    #fi

    currentTotal, currentBest = self.permutations[(nA,nU1,nU2)]

    minPermutationsNeeded = math.ceil((self.testCounter + 1) / 0.05)

    runningTotal = currentTotal
    runningPermutations = currentBest
    zscore = 0.0
    pvalue = 1.0
    nExceedences = bisect.bisect_left(runningPermutations, score)
    while runningTotal < minPermutationsNeeded:
        # If the number of exceedences is already greater than a certain number,
        # we can stop early, it is not significant
      if nExceedences > 10:
        break
      #fi

        # Otherwise, we can generate more permutations
      #newPermutations = self.parPool.map(self.generatePermutationScore,
      print("PERFORMING MORE PERMUTATIONS (%d/%d)" % (runningTotal, minPermutationsNeeded))
      newPermutations = [ self.generatePermutationScore(nA, nU1, nU2) for x in range(min(1000, minPermutationsNeeded - runningTotal)) ]
      runningTotal += 1000
      runningPermutations = sorted(runningPermutations + newPermutations)[:-1000]
      nExceedences = bisect.bisect_left(runningPermutations, score)
    #ewhile

      # Add the new permutations, if we made new ones
    if runningTotal > currentTotal:
      self.permutations[(nA, nU1, nU2)] = ( runningTotal, runningPermutations )
    #fi

      # Now that we have enough permutations, we can calculate a pvalue
    pvalue = 1 - float(runningTotal - nExceedences) / runningTotal
    zscore = stats.norm.ppf(1 - pvalue)


    return (zscore, pvalue)
  #edef

  ###############################################################################

  def generatePermutationScore(self, nA, nU1, nU2):
    import random
    A  = sum(random.sample(self.A, nA))
    U1 = sum(random.sample(self.U1, nU1))
    U2 = sum(random.sample(self.U2, nU2))
    return 2*A - (U1 + U2)
  #edef

  ###############################################################################

  def increaseTestCounter(self):
    self.testCounter += 1
    print(self.testCounter)
  #edef

#eclass
###############################################################################

class NodeTestInfo(object):

  #nGenes1 = None
  #nGenes2 = None
  #nA = None
  #nU1 = None
  #nU2 = None
  #cltTest = None
  #score = None
  #zscore = None
  #pvalue = None
  #lastNTests = None

  ###############################################################################

  def __init__(self, nGenes1, nGenes2, conservationScore, nA, nU1, nU2, score, cltTest, zscore, pvalue, lastNTests):
    self.nGenes1 = nGenes1
    self.nGenes2 = nGenes2
    self.conservationScore = conservationScore
    self.nA = nA
    self.nU1 = nU1
    self.nU2 = nU2
    self.score = score
    self.cltTest = cltTest
    self.zscore = zscore
    self.pvalue = pvalue
    self.lastNTests = lastNTests
  #edef

  ###############################################################################

  def isSignificant(self, alpha, nTests):
    if (self.pvalue * nTests) < alpha:
      return True
    else:
      return False
    #fi
   #edef

  ###############################################################################

  def setpvalue(self, pvalue):
    self.pvalue = pvalue
  #edef
 
  def setzscore(self, zscore):
    self.zscore = zscore
  #edef

  def setlastNTests(self, lastNTests):
    self.lastNTests = lastNTests
  #edef

  ###############################################################################     

#eclass

def writeClusterInfo(PS, clusters, outfile):
  print("writing now!")
  with open(outfile, "w") as ofd:
    ofd.write("#chrA\tchrB\tclusterID\tstartA\tendA\tstartB\tendB\tscore\tpvalue\tgenesA\tgenesB\n")
    for (chrpair, nodeID) in clusters:
      nodeInfo = PS.getNodeInfo(chrpair, nodeID)
      ((genome1Start, genome1End), (genome2Start, genome2End)) = PS.getNodeRegions(chrpair, nodeID)
      (genome1Genes, genome2Genes) = PS.getNodeGenes(chrpair, nodeID)
      ofd.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%s\n" % (chrpair[0], chrpair[1], nodeID, genome1Start, genome1End, genome2Start, genome2End, nodeInfo.score, nodeInfo.pvalue * PS.testCounter, ','.join(genome1Genes), ','.join(genome2Genes)))
    #efor
  #ewith
