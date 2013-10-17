# Need:
# cluster chr

import sys;

import bisect;
from scipy.cluster import hierarchy;
import numpy as np;
from ibidas import *;

import cluster_null as null;
reload(null);
import util as util;
reload(util);

###############################################################################

def prep_exon_list(E):
  L = zip(*E.Get(_.chrid, _.start)());

  chrs = {};
  for (chrid, start) in L:
    if chrid not in chrs:
      chrs[chrid] = [];
    #fi
    chrs[chrid].append(start);
  #efor

  for chrid in chrs.keys():
    chrs[chrid].sort();
  #efor

  return chrs;
#edef

###############################################################################

def count_exons_in_reg(chrs, chrid, start, end):
  L = chrs[chrid];
  i = bisect.bisect_right(L, start);
  j = bisect.bisect_left(L, end);

  return j - i + 1;
#edef

###############################################################################

def test_score(PR, k, reg):
  H, O, scores = PR.hit_windows[(k['id_a'], k['id_b'], k['WS'])];

  chrs_a = prep_exon_list(PR.org_exons[k['id_a']]);
  chrs_b = prep_exon_list(PR.org_exons[k['id_b']]);

  for (gen, r) in zip([chrs_a, chrs_b], reg[1]):
    print count_exons_in_reg(gen, r[1], r[2], r[3]);
  #efor
#edef

###############################################################################

def calc_distances(hits):
  HC = util.chr_pair_group(hits);
  D = {};

  nk = len(HC.keys());
  print "Calculating distances. This will take a while!";

  for (i, hc_k) in enumerate(HC.keys()):
    print "\r%d/%d" % (i+1, nk),
    sys.stdout.flush();
    hc  = HC[hc_k];
    if len(hc) < 2:
      continue;
    #fi
    cdm = np.array(condenseddm(hits, hc), dtype=np.dtype('u8'));
    if len(cdm) > 2:
      D[hc_k] = cdm;
    #fi
  #efor
  print "";

  return (HC, D);

#edef

###############################################################################

def calc_dendrograms(HC, D, linkage_type='single'):
  linkage_types = { 'single'   : hierarchy.single,
                    'complete' : hierarchy.complete,
                    'average'  : hierarchy.average,
                    'weighted' : hierarchy.weighted,
                    'centroid' : hierarchy.centroid,
                    'median'   : hierarchy.median,
                    'ward'     : hierarchy.ward };
  T = {};
  print "Calculating linkages. This will take a while!";
  nk = len(D.keys());
  for (i, dk) in enumerate(D.keys()):
    print "\r%d/%d" % (i+1, nk),
    sys.stdout.flush();
    L = linkage_types[linkage_type](D[dk]);
    T[dk] = construct_dendrogram(L, HC[dk]);
  #efor

  return T;

#edef

###############################################################################

def construct_dendrogram(Z, HC):
  """T = construct dendrogram(Z, HC)
     Z:  Linkage
     HC: List of hit IDS in this set
     H:  Sorted list of hits
   returns:
     * T: An array of binary tree nodes, the last node in the list is the root.
          Each node is of the form
          [ child1, child2, distance_between_children, hits_in_cluster ]
  """

  tree = [ dendrogram_node(None, None, 0, set([hc])) for hc in HC ];
  
  for z in Z:
    tree.append(dendrogram_node(int(z[0]), int(z[1]), z[2], tree[int(z[0])][3] | tree[int(z[1])][3]));
  #efor

  return tree;

#edef


###############################################################################

def score_dendrogram_node(T, n, hits, chrs_a, chrs_b):
  """score = score_dendrogram_node(T, r, H, O, chrs_a, chrs_b)
     Inputs:
       T       The linkage tree
       r:      The ID of the node in the linkage tree
       H:      The output from PR.windows
       O:      The output from PR.windows
       chrs_a: The output of prep_exon_list() for organism 0
       chrs_b: The output of prep_exon_list() for organism 1
     Outputs:
       score: A score for the cluster at this node in the dendrogram
  """

  I, J, dist, ids = T[n];
  hits = [ hits[i] for i in ids ];
  hit_ex_a = [ (h[2], h[3]) for h in hits ];
  hit_ex_b = [ (h[5], h[6]) for h in hits ];

  start_a = min([ h[7] for h in hits]);
  end_a   = max([ h[8] for h in hits]);

  start_b = min([ h[9]  for h in hits]);
  end_b   = max([ h[10] for h in hits]);

  exc_a = count_exons_in_reg(chrs_a, h[1], start_a, end_a);
  exc_b = count_exons_in_reg(chrs_b, h[4], start_b, end_b);

    # How many exons are "extra", or unaccounted for?
  dex_a = exc_a - len(set(hit_ex_a));
  dex_b = exc_b - len(set(hit_ex_b));

    # How many ends do we have in total?
  ends_a = len(hit_ex_a) + dex_a;
  ends_b = len(hit_ex_b) + dex_b;

  bitscore_sum = sum([h[13] for h in hits]);
  score = (2 * bitscore_sum) / (ends_a + ends_b);

  return (score, exc_a, exc_b);

#edef

###############################################################################

def calc_clusters_height(T, hits, chrs_a, chrs_b, H):

  C = [];

  for k in T.keys():
    t = T[k];

    c = cut_dendrogram_height(t, hits, chrs_a, chrs_b, H);

    C = C + c;
  #efor

  return C;

#edef

###############################################################################

def cut_dendrogram_height(T, hits, chrs_a, chrs_b, H=0):
  """C = cut_dendrogram(T, H, O, chrs_a, chrs_b, t)
   Inputs
       T       The linkage tree
       hits:   The indexed hits
       chrs_a: The output of prep_exon_list() for organism 0
       chrs_b: The output of prep_exon_list() for organism 1
       D:      The depth to cut at.
     Outputs: 
       C: A list of clusters of the form (L, score) where L is a list of hits in the cluster, and score is the score associated with that cluster
  """

  istack = [ -1 ];
  C = [];
  D = {};

  while len(istack) > 0:

    index = istack.pop();

    I, J, dist, ids     = T[index];
    score, exc_a, exc_b = score_dendrogram_node(T, index, hits, chrs_a, chrs_b);

    if dist <= H:
      cdesc = clust_description(hits, ids, score);
      C.append(cdesc);
    else:
      if not(I is None):
        istack.append(I);
      #fi
      if not(J is None):
        istack.append(J);
      #fi
    #fi
  #efor

  return C;
#edef

###############################################################################

def calc_clusters(T, hits, chrs_a, chrs_b, t=None):

  C = [];

  for k in T.keys():
    t = T[k];

    c = cut_dendrogram(t, hits, chrs_a, chrs_b);

    C = C + c;
  #efor

  return C;

#edef

###############################################################################    

def cut_dendrogram(T, hits, chrs_a, chrs_b, t=0.05):
  """C = cut_dendrogram(T, H, O, chrs_a, chrs_b, t)
   Inputs
       T       The linkage tree
       hits:   The indexed hits
       chrs_a: The output of prep_exon_list() for organism 0
       chrs_b: The output of prep_exon_list() for organism 1
     Outputs: 
       C: A list of clusters of the form (L, score) where L is a list of hits in the cluster, and score is the score associated with that cluster
  """

  istack = [ -1 ];
  C = [];
  D = {};

  dsize = lambda d: sum([ len(d[k]) for k in d])

  nd = null.cluster_null_clt([ h[13] for h in hits], dsize(chrs_a), dsize(chrs_b));
                              

  while len(istack) > 0:

    index = istack.pop();

    I, J, dist, ids     = T[index];
    score, exc_a, exc_b = score_dendrogram_node(T, index, hits, chrs_a, chrs_b);

    p = nd.pvalue(exc_a, exc_b, score);
    print '%2.10f %10d %10d %2.10f' % (p, exc_a, exc_b, score);
    if p < 0.05:
      cdesc = clust_description(hits, ids, score, p);
      C.append(cdesc);
    else:
      if not(I is None):
        istack.append(I);
      #fi
      if not(J is None):
        istack.append(J);
      #fi
    #fi
  #efor

  return C;
#edef
  

###############################################################################

def dendrogram_node(i, j, dist, ids):
  return (i, j, dist, ids);
#edef

###############################################################################

def distance(hits, i, j):
  """d = distance(hits, i, j):
     hits: The output list of hit_index
     i: The ID of hit i
     j: The ID of hit j

                |-a-|
     ---=========---============---
        \\\\\\\\\\  ||||||||||||
     ----=========--============---
                 |-b|

      d = a + b

     Outputs:
       d: The distance between hit i and hit j
  """
  h1 = hits[i];
  h2 = hits[j];

    # Different chromosomes
  if (h1[1] != h2[1]) or (h1[4] != h2[4]):
    return float("inf");
  #fi

    # Same exons
  if (h1[2] == h2[2]) and (h1[5] == h2[5]) and (h1[3] == h2[3]) and (h1[6] == h2[6]):
    return 0;
  #fi

    # Regions of hit
  h1_a = (h1[7], h1[8]);
  h1_b = (h1[9], h1[10]);
  h2_a = (h2[7], h2[8]);
  h2_b = (h2[9], h2[10]);
  ov1 = util.overlap(h1_a, h2_a);
  ov2 = util.overlap(h1_b, h2_b);

    # If they all overlap, distance is 0.
  if ov1 > 0 and ov2 > 0:
    return 0;
  #fi

  return -(min(0, ov1) + min(0, ov2));
#edef

###############################################################################

def condenseddm(hits, L):
  """cdm = condenseddm(H, O, L)
     H: Output list from sort_org
     O: Output list from sort_org
     L: A list of hits

     Outputs:
       cdm:  A condensed distance matrix (Upper triangle of distance matrix as array
  """

  nl = len(L);

  cdm = [];

  for i in xrange(nl-1):
    for j in xrange(i+1, nl):
      cdm.append(distance(hits, L[i], L[j]));
    #efor
  #efor

  return cdm;
#edef

###############################################################################

def clust_description(hits, C, score, p):
  """cd = clust_description(H, O, scores, C):
     H:     Output list of sort_org
     O:     Output list of sort_org
     score: Scores per hit
     p:     The pvalue of this cluster
     C:     Hits in cluster

     Outputs:
       cd: A description of the cluster.
  """

  hits   = [ hits[int(i)] for i in C ];
  n_hits = len(hits);

  prots_a = set([ h[2] for h in hits ]);
  prots_b = set([ h[5] for h in hits ]);

  a_chr   = hits[0][1];
  a_start = min([ h[7] for h in hits ]);
  a_end   = max([ h[8] for h in hits ]);

  b_chr   = hits[0][4];
  b_start = min([ h[9]  for h in hits ]);
  b_end   = max([ h[10] for h in hits ]);

  return (a_chr, a_start, a_end, b_chr, b_start, b_end, n_hits, score, p, prots_a, prots_b, C);
#edef

###############################################################################
