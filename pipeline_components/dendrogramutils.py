
from collections import namedtuple
import csv

import utils as utils

nodeType = namedtuple("dendrogramNode", "child1 child2 leaves")

###############################################################################

def dendrogramNode(child1, child2, children):
  return nodeType( child1, child2, set(children))
#edef

###############################################################################

def prepareDendrogram(oD):

  nObjects = len(oD) + 1

  D = [ dendrogramNode(-1, -1, set([nodeI])) for nodeI in range(nObjects) ]

  for node in oD:
    D.append(dendrogramNode( int(node[0]), int(node[1]), D[int(node[0])].leaves | D[int(node[1])].leaves ) )
  #efor
  return D
#edef
  

###############################################################################

def readDendrogram(dFile):
  D = []
  with open(dFile, "r") as ifd:
    reader = csv.reader(ifd, delimiter='\t')
    for row in reader:
      if row[0][0] == '#':
        continue;
      #fi
      D.append( dendrogramNode(int(row[0]), int(row[1]), set([ int(c) for c in row[2].split(",")])) )
    #efor
  #ewith
  return D
#edef

###############################################################################

def writeDendrogram(D, dFile):
  with open(dFile, "w") as ofd:
    ofd.write("#child1\tchild2\tleaves\n")
    for node in D:
      ofd.write("%d\t%d\t%s\n" % (node.child1, node.child2, ','.join([ str(x) for x in node.leaves ])))
    #efor
  #ewith
#edef
