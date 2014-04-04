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

from ibidas.utils.util import debug_here;

###############################################################################

def null_dist(dist=null.cluster_null_score_strict_smart, **kwargs):
    return dist(**kwargs);
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

def score_dendrogram_node(T, n, hits, chrs_a, chrs_b, hitA, hitB):
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
  H = [ hits[i] for i in ids ];

  return score_hits(H, hits, chrs_a, chrs_b, hitA, hitB);
#edef

###############################################################################

def score_hits(H, hits, chrs_a, chrs_b, hitA, hitB):
  """Given a set of hits, calculate a cluster score."""

  print "----START SCORE"

  hit_ex_a = [ (h[2], h[3]) for h in H ];
  hit_ex_b = [ (h[5], h[6]) for h in H ];

  start_a = min([ h[7] for h in H]);
  end_a   = max([ h[8] for h in H]);

  start_b = min([ h[9]  for h in H]);
  end_b   = max([ h[10] for h in H]);

  ex_a = util.exons_in_reg(chrs_a, h[1], start_a, end_a);
  ex_b = util.exons_in_reg(chrs_b, h[4], start_b, end_b);

  UE_a = set(ex_a) - set(hit_ex_a);
  UE_b = set(ex_b) - set(hit_ex_b);

  nz_UE_a = [ k for k in UE_a if k in hitA ];
  nz_UE_b = [ k for k in UE_b if k in hitB ];

  SE_a = sum([ hits[hitA[k]][12] for k in nz_UE_a ]);
  SE_b = sum([ hits[hitB[k]][12] for k in nz_UE_b ]);

  score_sum = sum([h[12] for h in H]);
  score = 2*score_sum - SE_a - SE_b;

  return (score, ex_a, ex_b, UE_a, UE_b, nz_UE_a, nz_UE_b);

#edef

###############################################################################

def calc_clusters(T, hits, chrs_a, chrs_b, cut, alpha, dist, ngenes_threshold, conservation_ratio):
  """Detect significant clusters based on the dynamic cutting algorithm."""


  C = [];

    # Precompute max hit scores per exon
  hitA = {};
  hitB = {};
  v     = [ (2,3,hitA), (5,6,hitB) ];

  for (i, h) in enumerate(hits):
    for (j, k, H) in v:
      k = (h[j], h[k]);
      if not(k in H):
        H[k] = i;
      elif h[12] > hits[H[k]][12]:
         H[k] = i;
      #fi
    #efor
  #efor

  tests = 0;
            # Chromosome 1, chromosome2, index, score, pvalue)
  ntrees       = len(T);
  start_clusts = [ (k[0], k[1], -1, None, None) for k in T.keys() ];
  #start_clusts = [ (1,1,-1,None,None)];
  sig_clusts   = [];
  C            = [];

  nd = null_dist(dist=dist,
                 alpha = alpha,
                 tests = len(start_clusts),
                 scores=[ h[12] for h in hits],
                 hits=hits, hitA=hitA, hitB=hitB,
                 storage='null_dist_store.dat');

  # Initial loop through all trees individually
  for (i, sc) in enumerate(start_clusts):
    print "\nTree %d/%d" % (i+1, len(start_clusts));
    nd.update_mperm(len(start_clusts));
    sig_clusts = sig_clusts + complete_set(T, [ sc ], hits, chrs_a, chrs_b, hitA, hitB, nd, cut, ngenes_threshold, conservation_ratio);
    tests = tests + nd.tests - (ntrees - 1);
    print tests;
  #efor

  nd.update_mperm(tests);
  print "Tests done: %d\nSignificant Clusters: %d" % (nd.tests, len(sig_clusts));

  # Recalculate p-values for clusters deemed significant
  sig_clusts = complete_set(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, cut, ngenes_threshold, conservation_ratio);

  for (c1, c2, index, score, pvalue) in sig_clusts:
    I, J, dist, ids     = T[(c1,c2)][index];
    cdesc = clust_description(hits, ids, score, pvalue);
    C.append(cdesc);
  #efor

  return C;

#edef

###############################################################################

def complete_set(T, start_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, cut, ngenes_threshold, conservation_ratio):
  """Perform one loop of the dynamic correction algorithm."""

  tests = nd.tests;

  sig_clusts = start_clusts

  while True:
    if cut == 'simple':
      sig_clusts = cut_dendrogram_tests_simple(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd);
    elif cut == 'greater':
      sig_clusts = cut_dendrogram_tests_greater(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd);
    elif cut == 'deeper_greater':
      sig_clusts = cut_dendrogram_tests_deeper_greater(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold, conservation_ratio);
    elif cut == 'deeper':
      sig_clusts = cut_dendrogram_tests_deeper_greater(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold, conservation_ratio=0.0)
    elif cut == 'max':
      return cut_dendrogram_tests_max(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold);
    elif cut == 'peaks':
      return cut_dendrogram_tests_peaks(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold);
    else:
      print "Unknown cut type: %s" % (cut);
      return [];
    
    if nd.tests == tests:
      break;
    else:
      print "We need to loop through again!"
      tests = nd.tests;
    #fi
  #ewhile

  return sig_clusts;
#edef

###############################################################################

def cut_dendrogram_tests_simple(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold=2):
  """ Cut as soon as we find a significant cluster"""

  iqueue = sig_clusts;
  SC     = [];

  i = 0;

  while len(iqueue) > 0:
    print "iqueue: %d, sigclusts: %d, done: %d" % (len(iqueue), len(SC), i+1);
    (c1, c2, index, C_info, old_pvalue) = iqueue.pop();
    t               = T[(c1,c2)];
    I, J, dist, ids = t[index];
    i = i + 1;

    if C_info == None:
      C_info = score_dendrogram_node(t, index, hits, chrs_a, chrs_b, hitA, hitB);
      C_info = (C_info[0], len(ids), len(C_info[5]), len(C_info[6]));
    #fi
    score, n, nz_nue_a, nz_nue_b = C_info;

    pvalue = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);

    if pvalue < (nd.alpha / nd.tests):
      SC.append((c1, c2, index, C_info, pvalue));
      print "SIGNIFICANT!!!!", SC[-1];
      continue;
    else:
      additional_tests = 0;
      I_ngenes_a = len(set([hits[h][2] for h in t[I][3]]));
      I_ngenes_b = len(set([hits[h][5] for h in t[I][3]]));
      suff_I = (I_ngenes_a > ngenes_threshold) or (I_ngenes_b > ngenes_threshold);

      J_ngenes_a = len(set([hits[h][2] for h in t[J][3]]));
      J_ngenes_b = len(set([hits[h][5] for h in t[J][3]]));
      suff_J = (J_ngenes_a > ngenes_threshold) or (J_ngenes_b > ngenes_threshold);

      if I is not None and suff_I:
        iqueue.insert(0, (c1, c2, I, None, None));
        additional_tests = additional_tests + 1;
      #fi
      if J is not None and suff_J:
        iqueue.insert(0, (c1, c2, J, None, None));
        additional_tests = additional_tests + 1;
      #fi
      nd.update_mperm(nd.tests + additional_tests);
    #fi
  #ewhile

  print "TD: %d, SC: %d TTD: %d\r" %(nd.tests, len(SC), len(iqueue))
  return SC;

#edef


###############################################################################

def cut_dendrogram_tests_greater(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold=2, conservation_ratio=1):
  """Cut as soon as we find a significant cluster that fulfills the conservation_ratio"""

  iqueue = sig_clusts;
  SC     = [];

  i = 0;

  while len(iqueue) > 0:
    print "iqueue: %d, sigclusts: %d, done: %d" % (len(iqueue), len(SC), i+1);
    (c1, c2, index, C_info, old_pvalue) = iqueue.pop();
    t               = T[(c1,c2)];
    I, J, dist, ids = t[index];
    i = i + 1;

    if C_info == None:
      C_info = score_dendrogram_node(t, index, hits, chrs_a, chrs_b, hitA, hitB);
      C_info = (C_info[0], len(ids), len(C_info[5]), len(C_info[6]));
    #fi
    score, n, nz_nue_a, nz_nue_b = C_info;

      # Use only clusters which have more hits than unaccounted hits
    
    if float(n)/float(nz_nue_a + nz_nue_b) < conservation_ratio:
      pvalue = 1.0;
    else:
      pvalue = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);
    #fi

    if pvalue < (nd.alpha / nd.tests):
      SC.append((c1, c2, index, C_info, pvalue));
      print "SIGNIFICANT!!!!", SC[-1];
      continue;
    else:
      additional_tests = 0;
      I_ngenes_a = len(set([hits[h][2] for h in t[I][3]]));
      I_ngenes_b = len(set([hits[h][5] for h in t[I][3]]));
      suff_I = (I_ngenes_a > ngenes_threshold) or (I_ngenes_b > ngenes_threshold);

      J_ngenes_a = len(set([hits[h][2] for h in t[J][3]]));
      J_ngenes_b = len(set([hits[h][5] for h in t[J][3]]));
      suff_J = (J_ngenes_a > ngenes_threshold) or (J_ngenes_b > ngenes_threshold);

      if I is not None and suff_I:
        iqueue.insert(0, (c1, c2, I, None, None));
        additional_tests = additional_tests + 1;
      #fi
      if J is not None and suff_J:
        iqueue.insert(0, (c1, c2, J, None, None));
        additional_tests = additional_tests + 1;
      #fi
      nd.update_mperm(nd.tests + additional_tests);
    #fi
  #ewhile

  print "TD: %d, SC: %d TTD: %d\r" %(nd.tests, len(SC), len(iqueue))
  return SC;

#edef


###############################################################################

def cut_dendrogram_tests_max(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold=2):
  """Descend into the tree until we find the maximum node in the path. Return all maximum nodes in the branch"""

  if len(sig_clusts) == 1:
    istack = [ (a,b,c,None,(0.0,1.0),0.0,0) for (a,b,c,d,e) in sig_clusts ];
    SC     = [];

    i_count = 0;
    while len(istack) > 0:
      (c1, c2, index, C_info, (zscore,pvalue), max_zscore, have_sc) = istack.pop();
      t               = T[(c1,c2)];
      I, J, dist, ids = t[index];

      ngenes_a = len(set([hits[h][2] for h in t[index][3]]));
      ngenes_b = len(set([hits[h][5] for h in t[index][3]]));
      suff     = (ngenes_a > ngenes_threshold) or (ngenes_b > ngenes_threshold);
   
      if not(suff):
        continue;
      #fi

      if C_info == None:
        C_info_r = score_dendrogram_node(t, index, hits, chrs_a, chrs_b, hitA, hitB);
        C_info   = (C_info_r[0], len(ids), len(C_info_r[5]), len(C_info_r[6]));

        nd.update_mperm(nd.tests + 1);
        score, n, nz_nue_a, nz_nue_b = C_info;
        (zscore, pvalue)             = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);

        max_zscore = max(max_zscore, zscore);
        istack.append((c1, c2, index, C_info, (zscore,pvalue), max_zscore, len(SC)));
        istack.append((c1, c2, I,     None,   (0.0,1.0),       max_zscore, len(SC)));
        istack.append((c1, c2, J,     None,   (0.0,1.0),       max_zscore, len(SC)));

        print len(istack), len(t), i_count;
        i_count = i_count + 1;
      else:
        print "WE HAVE RETURNED!";
        # We have returned!!
        if len(SC) == have_sc and zscore >= max_zscore:
          # We did not find any significant clusters beneath this node, and this IS the significant node!
          print "WE FOUND MAX SIG!";
          SC.append((c1,c2,index,C_info,pvalue));
        #fi
    #ewhile

  else:
    SC = sig_clusts
  #fi

  SCF = [];
  for sc in SC:
    (c1,c2,index,C_info,pvalue)  = sc;
    score, n, nz_nue_a, nz_nue_b = C_info

    # Recalculate the pvalue here
    (zscore,pvalue) = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);
    if pvalue < (nd.alpha / float(len(SC))):
      SCF.append(sc);
    #fi
  #efor

  print SCF;

  return SCF;

#edef

###############################################################################

def cut_dendrogram_tests_peaks(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold=2):
  """Return all clusters which are more significant than their children and their parents."""

  if len(sig_clusts) == 1:
    istack = [ (a,b,c,d,e,False) for (a,b,c,d,e) in sig_clusts ];
    SC     = [];
    
    i_count = 0;
    while len(istack) > 0:
      clust = istack.pop();
      (c1, c2, index, C_info, measure, parent_zscore) = clust;
      t               = T[(c1,c2)];
      I, J, dist, ids = t[index];
      
      ngenes_a = len(set([hits[h][2] for h in t[index][3]]));
      ngenes_b = len(set([hits[h][5] for h in t[index][3]]));
      suff     = (ngenes_a > ngenes_threshold) or (ngenes_b > ngenes_threshold);
      
      if not(suff):
        continue;
      #fi

      if C_info == None:
        C_info_r = score_dendrogram_node(t, index, hits, chrs_a, chrs_b, hitA, hitB);
        C_info   = (C_info_r[0], len(ids), len(C_info_r[5]), len(C_info_r[6]));

        score, n, nz_nue_a, nz_nue_b = C_info;
        (zscore, pvalue)             = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);
      else:
        (zscore, pvalue) = measure;
      #fi
      
      (I_zscore, I_pvalue) = (0.0, 1.0);
      (J_zscore, J_pvalue) = (0.0, 1.0);

      if I != None:
        I_info_r             = score_dendrogram_node(t, I, hits, chrs_a, chrs_b, hitA, hitB);
        I_info               = (I_info_r[0], len(t[I][3]), len(I_info_r[5]), len(I_info_r[6]));
        (I_zscore, I_pvalue) = nd.pvalue(score=I_info[0], n=I_info[1], nue_a=I_info[2], nue_b=I_info[3]);
        istack.append((c1, c2, I, I_info, (I_zscore, I_pvalue), zscore));
      #fi

      if J != None:
        J_info_r             = score_dendrogram_node(t, J, hits, chrs_a, chrs_b, hitA, hitB);
        J_info               = (J_info_r[0], len(t[J][3]), len(J_info_r[5]), len(J_info_r[6]));
        (J_zscore, J_pvalue) = nd.pvalue(score=J_info[0], n=J_info[1], nue_a=J_info[2], nue_b=J_info[3]);
        istack.append((c1, c2, J, J_info, (J_zscore, J_pvalue), zscore));
      #fi

      if ((zscore > I_zscore) or (zscore > J_zscore)) and zscore > parent_zscore:
        print "WE FOUND A PEAK!";
        SC.append((c1, c2, index, C_info, (zscore,pvalue)));
      #fi

      

      print len(istack), len(t), i_count;
      i_count = i_count + 1;
    #ewhile

  else:
    SC = sig_clusts
  #fi

  SCF = [];
  for sc in SC:
    (c1,c2,index,C_info,pvalue)  = sc;
    score, n, nz_nue_a, nz_nue_b = C_info
    
    # Recalculate the pvalue here
    (zscore,pvalue) = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);
    if pvalue < (nd.alpha / float(len(SC))):
      SCF.append(sc);
    #fi
  #efor

  print SCF;

  return SCF;

#edef

###############################################################################

def cut_dendrogram_tests_deeper(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold=2):
  """Descend to the child nodes if the parent is less significant than the parent."""

  iqueue = sig_clusts;
  SC     = [];

  i = 0;

  while len(iqueue) > 0:
    print "iqueue: %d, sigclusts: %d, done: %d" % (len(iqueue), len(SC), i+1);
    (c1, c2, index, C_info, old_pvalue) = iqueue.pop();
    t               = T[(c1,c2)];
    I, J, dist, ids = t[index];
    i = i + 1;
    deeper = 0; # We have not found any /more/ significant clusters... yet

      # Only calculate the score if we don't already have it...
    if C_info == None:
      C_info_r = score_dendrogram_node(t, index, hits, chrs_a, chrs_b, hitA, hitB);
      C_info = (C_info_r[0], len(ids), len(C_info_r[5]), len(C_info_r[6]));
    #fi
    score, n, nz_nue_a, nz_nue_b = C_info;

      # Check if children would be eligible for descending into.
    suff_I = I is not None;
    if suff_I:
      I_ngenes_a = len(set([hits[h][2] for h in t[I][3]]));
      I_ngenes_b = len(set([hits[h][5] for h in t[I][3]]));
      suff_I = (I_ngenes_a > ngenes_threshold) or (I_ngenes_b > ngenes_threshold);
    #fi

    suff_J = J is not None;
    if suff_J:
      J_ngenes_a = len(set([hits[h][2] for h in t[J][3]]));
      J_ngenes_b = len(set([hits[h][5] for h in t[J][3]]));
      suff_J = (J_ngenes_a > ngenes_threshold) or (J_ngenes_b > ngenes_threshold);
    #fi

      # Get pvalue for cluster
    (zscore, pvalue) = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);

      # OK, we have a significant cluster. Now check if it's children are /more/ significant...
    if pvalue < (nd.alpha / nd.tests):
      (I_zscore, I_pvalue) = (0.0, 1.0);
      (J_zscore, J_pvalue) = (0.0, 1.0);

        # Calculate pvalues for child I
      if suff_I:
        I_info_r = score_dendrogram_node(t, I, hits, chrs_a, chrs_b, hitA, hitB);
        I_info   = (I_info_r[0], len(t[I][3]), len(I_info_r[5]), len(I_info_r[6]));
        (I_zscore, I_pvalue) = nd.pvalue(score=I_info[0], n=I_info[1], nue_a=I_info[2], nue_b=I_info[3]);
      #fi

        # calculate pvalues for child J
      if suff_J:
        J_info_r = score_dendrogram_node(t, J, hits, chrs_a, chrs_b, hitA, hitB);
        J_info   = (J_info_r[0], len(t[J][3]), len(J_info_r[5]), len(J_info_r[6]));
        (J_zscore, J_pvalue) = nd.pvalue(score=J_info[0], n=J_info[1], nue_a=J_info[2], nue_b=J_info[3]);
      #fi

      if I_zscore > zscore:
        iqueue.append((c1, c2, I, I_info, I_pvalue));
        print "DESCENDING LEFT: MORE SIG";
        deeper = deeper + 1;
      #fi

      if J_zscore > zscore:
        iqueue.append((c1, c2, J, J_info, J_pvalue));
        print "DESCENDING RIGHT: MORE SIG";
        deeper = deeper + 1;
      #fi

      if deeper == 0:
        SC.append((c1, c2, index, C_info, pvalue));
        print "SIGNIFICANT!!!!", SC[-1];
      #fi
    else:

      if suff_I:
        iqueue.append((c1, c2, I, None, None));
        print "DESCENDING LEFT";
        deeper = deeper + 1;
      #fi
      if suff_J:
        iqueue.append((c1, c2, J, None, None));
        print "DESCENDING RIGHT";
        deeper = deeper + 1;
      #fi
    #fi

    nd.update_mperm(nd.tests + deeper);
  #ewhile

  print "TD: %d, SC: %d TTD: %d\r" %(nd.tests, len(SC), len(iqueue))
  return SC;

#edef

###############################################################################


def cut_dendrogram_tests_deeper_greater(T, sig_clusts, hits, chrs_a, chrs_b, hitA, hitB, nd, ngenes_threshold, conservation_ratio):
  """Descend into the child node if the parent is less significant or doesn't satisfy the conservation ratio."""


  iqueue = sig_clusts;
  SC     = [];

  i = 0;

  while len(iqueue) > 0:
    print "iqueue: %d, sigclusts: %d, done: %d" % (len(iqueue), len(SC), i+1);
    (c1, c2, index, C_info, old_pvalue) = iqueue.pop();
    t               = T[(c1,c2)];
    I, J, dist, ids = t[index];
    i = i + 1;
    deeper = 0; # We have not found any /more/ significant clusters... yet

      # Only calculate the score if we don't already have it...
    if C_info == None:
      C_info_r = score_dendrogram_node(t, index, hits, chrs_a, chrs_b, hitA, hitB);
      C_info = (C_info_r[0], len(ids), len(C_info_r[5]), len(C_info_r[6]));
    #fi
    score, n, nz_nue_a, nz_nue_b = C_info;

      # Check if children would be eligible for descending into.
    suff_I = I is not None;
    if suff_I:
      I_ngenes_a = len(set([hits[h][2] for h in t[I][3]]));
      I_ngenes_b = len(set([hits[h][5] for h in t[I][3]]));

      I_info_r = score_dendrogram_node(t, I, hits, chrs_a, chrs_b, hitA, hitB);
      I_info   = (I_info_r[0], len(t[I][3]), len(I_info_r[5]), len(I_info_r[6]))
      suff_I   = ((I_ngenes_a > ngenes_threshold) or (I_ngenes_b > ngenes_threshold)) #and (float(I_info[1] + 1) / float(I_info[2] + I_info[3] + 1)) >= conservation_ratio;
    #fi

    suff_J = J is not None;
    if suff_J:
      J_ngenes_a = len(set([hits[h][2] for h in t[J][3]]));
      J_ngenes_b = len(set([hits[h][5] for h in t[J][3]]));

      J_info_r = score_dendrogram_node(t, J, hits, chrs_a, chrs_b, hitA, hitB);
      J_info   = (J_info_r[0], len(t[J][3]), len(J_info_r[5]), len(J_info_r[6]));
      suff_J   = ((J_ngenes_a > ngenes_threshold) or (J_ngenes_b > ngenes_threshold)) #and (float(J_info[1] + 1) / float(J_info[2] + J_info[3] + 1)) >= conservation_ratio;
    #fi

      # Get pvalue for cluster
    ratio = float(n+1) / float(nz_nue_a + nz_nue_b + 1);
    if ratio >= conservation_ratio:
      (zscore, pvalue) = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);
    else:
      (zscore, pvalue) = (0.0, 1.0);
    #fi

      # OK, we have a significant cluster. Now check if it's children are /more/ significant...
    if pvalue < (nd.alpha / nd.tests):
      (I_zscore, I_pvalue) = (0.0, 1.0);
      (J_zscore, J_pvalue) = (0.0, 1.0);

        # Calculate pvalues for child I
      if suff_I and (float(I_info[1] + 1) / float(I_info[2] + I_info[3] + 1)) >= conservation_ratio:
        (I_zscore, I_pvalue) = nd.pvalue(score=I_info[0], n=I_info[1], nue_a=I_info[2], nue_b=I_info[3]);
      #fi

        # calculate pvalues for child J
      if suff_J and (float(J_info[1] + 1) / float(J_info[2] + J_info[3] + 1)) >= conservation_ratio:
        (I_zscore, I_pvalue) = nd.pvalue(score=J_info[0], n=J_info[1], nue_a=J_info[2], nue_b=J_info[3]);
      #fi

        # Due to fp limitations, we compare zscores instead of pvalues here
      if (I_zscore > zscore) or (J_zscore > zscore):
        iqueue.append((c1, c2, I, I_info, I_pvalue));
        iqueue.append((c1, c2, J, J_info, J_pvalue));
        print "DESCENDING: MORE SIG";
        deeper = deeper + 2;
      else:
        SC.append((c1, c2, index, C_info, pvalue));
        print "SIGNIFICANT!!!!", SC[-1];
      #fi
    else:

      if suff_I:
        iqueue.append((c1, c2, I, None, None));
        print "DESCENDING LEFT";
        deeper = deeper + 1;
      #fi
      if suff_J:
        iqueue.append((c1, c2, J, None, None));
        print "DESCENDING RIGHT";
        deeper = deeper + 1;
      #fi
    #fi

    nd.update_mperm(nd.tests + deeper);
  #ewhile

  print "TD: %d, SC: %d TTD: %d\r" %(nd.tests, len(SC), len(iqueue))
  return SC;

#edef


###############################################################################    

def pi(T, index, hits, chrs_a, chrs_b, hitA, hitB, nd):
  if (index is None):
    return (None, None, None, None, None, 1.0);
  #fi

  I, J, dist, ids               = T[index];
  score, ex_a, ex_b, ue_a, ue_b, nz_ue_a. nz_ue_b = score_dendrogram_node(T, index, hits, chrs_a, chrs_b, hitA, hitB);

  p = nd.pvalue(score=score, n=len(ids), nue_a=len(ue_a), nue_b=len(ue_b));
    

  return (score, ex_a, ex_b, ue_a, ue_b, p);
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
       \\\\\\\\\\   ||||||||||||
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

