import sys;
sys.path.append('./utils');

import proteny as ps;
from utils import cluster;
from utils import cluster_null as null;
from utils import util;
reload(cluster);
reload(null);

import matplotlib.pyplot as plt;
import pylab;

import math;
import util;

import bisect;

from scipy.cluster import hierarchy;
import numpy as np;

###############################################################################

def overlap_clusters(C):

  overlaps = [];

  for i in xrange(len(C)-1):
    ci  = C[i];
    cio = [];
    for j in xrange(i+1, len(C)):
      cj = C[j];
      rlen = float(ci[2] - ci[1] + 1) / float(cj[2] - cj[1] + 1);
      if ((ci[0] == cj[0]) and        
          util.overlap( (ci[1], ci[2]), (cj[1], cj[2])) > 0 and
          (rlen > 0.8 and rlen < 1.25)):
        cio.append(j);
      #fi
    #efor
    if len(cio) > 0:
      overlaps.append(cio + [i]);
    #fi
  #efor

  return overlaps;

#edef

###############################################################################

def overlap_clusters2(C):
  overlaps = [];

  for i in xrange(len(C)):
    print "\r%d / %d" % (i, len(C)),
    sys.stdout.flush()
    K = [ C[i] ];
    for j in xrange(i+1, len(C)):
      #fi
      cj = C[j];
      for k in xrange(len(K)):
        ck = K[k];
        len1 = float(cj[2] - cj[1] + 1) / float(ck[2] - ck[1] + 1);
        len2 = float(cj[5] - cj[4] + 1) / float(ck[5] - ck[4] + 1); 

          # If one of the two regions
          # * overlap, and
          # * have a reasonable similar size,
          # add it to the overlap region
        if (((cj[0] == ck[0]) and
             util.overlap( (cj[1], cj[2]), (ck[1], ck[2])) > 0 and
             (len1 > 0.3 and len1 < 3))
            or
            ((cj[3] == ck[3]) and
             util.overlap( (cj[4], cj[5]), (ck[4], ck[5])) > 0 and
             (len2 > 0.3 and len2 < 3))):
          K.append(cj);
          break;
        #fi
      #efor
    #efor
    if len(K) > 1:
      overlaps.append(K);
    #fi
  #efor

    # Overlapping regions
  OR = [];

    # Within an overlap, find unique regions; remove duplicate regions
  for k in overlaps:
    R1 = [ (c[0], c[1], c[2], 'id_a') for c in k ];
    R2 = [ (c[3], c[4], c[5], 'id_b') for c in k ];

      # List of unique regions
    UR = [];

    for R in [ R1, R2 ]:
      for i in xrange(len(R)-1):
        mr = R[i]; # maximal region
        if mr == None:
          continue;
        for j in xrange(i+1, len(R)):
          rj = R[j];
          if rj == None:
            continue;
            # If the two regions overlap, expand the maximal region
          if (mr[0] == rj[0] and
              util.overlap( (mr[1], mr[2]), (rj[1], rj[2])) > 0 ):
            mr = (mr[0], min(mr[1], rj[1]), max(mr[2], rj[2]), mr[3]);
            R[j] = None;
          #fi
        #efor
        UR.append(mr);
      #efor
    #efor
    OR.append(UR);
  #efor

    # Remove duplicates among overlaps (check for subsets)
  for i in xrange(len(OR)-1):
    for j in xrange(i+1, len(OR)):
      if sum([ 1 for k in OR[j] if (k in OR[i]) ]) >= 0.8 * len(OR):
        if len(OR[i]) > len(OR[j]):
          OR[j] = [];
        else:
          OR[i] = [];
      #fi
    #efor
  #efor
         

  return [ x for x in OR if len(x) > 0 ];

#edef

###############################################################################

def hit_cluster_overlap2(PR, k):
  """ Identify clusters which are overlapped"""

  if ((k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd']) not in PR.hit_clusters):
    print "You must run hit_cluster() first!";
    return None;
  #fi

  C = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd'])];


  OR = overlap_clusters2(C);
  R  = [];

  for (j, UR) in enumerate(OR):
    regs = [];
    for r in UR:
      regs.append((PR.org_names[k[r[3]]], r[0], r[1], r[2]));
    #efor
    R.append( ( str(j), regs));
  #efor

  return sorted(R,key=lambda x:len(x[1]));

 #efor

###############################################################################

def hit_cluster_overlap(PR, k):
  """ Identify clusters which are overlapped"""

  if ((k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd']) not in PR.hit_clusters):
    print "You must run hit_cluster() first!";
    return None;
  #fi

  C = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd'])];

  overlaps = overlap_clusters(C);
  R        = [];

  for (j, cio) in enumerate(overlaps):
    areg_chr   = C[cio[0]][0];
    areg_start = min([ C[i][1] for i in cio ]);
    areg_end   = max([ C[i][2] for i in cio ]);
    areg       = ( PR.org_names[k['id_a']], areg_chr, areg_start, areg_end );

    bregs      = [ ( PR.org_names[k['id_b']], C[i][3], C[i][4], C[i][5] ) for i in cio ];
    regs       = [ areg ] + bregs;
    R.append( ( str(j), regs));
  #efor

  return sorted(R,key=lambda x:x[1]);

 #efor

###############################################################################

def hit_clusters_containing(PR, k, a_prots, b_prots):


  C  = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd'])];
  CA = [];
  CB = [];

  for ap in a_prots:
    CA.append([ i for (i,c) in enumerate(C) if ap in c[9] ]);
  #efor
  for bp in b_prots:
    CB.append([ i for (i,c) in enumerate(C) if bp in c[10] ]);
  #efor

  return (CA, CB);
#edef

###############################################################################

def hit_clust_2_reg(PR, k, ci, name=None):
  C  = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd'])][ci];

  #('bri1',   [ ('schco2', 'scaffold_1',  5346259, 5348750), ('agabi', 'scaffold_1',  1865954, 1869064) ] ),

  reg = (name if not(name == None) else str(ci),
         [ ( PR.org_names[k['id_a']], C[0], C[1], C[2] ),
           ( PR.org_names[k['id_b']], C[3], C[4], C[5] ) ] );

  return reg;
#edef

###############################################################################

def delndups(PR, k):

  """ Return a list of clusters which have an obvious deletion or duplication """

  dups = [ i for (i,C) in enumerate(PR.hit_clusters[k])  if not(len(C[9]) == len(C[10])) ];

  return dups;

#edef

###############################################################################

def cmp_clt_perm(PR, k):

  hits  = PR.hits[(k['id_a'], k['id_b'])];
  T     = PR.hit_dendrograms[(k['id_a'], k['id_b'], k['linkage_type'])];

  chrs_a = PR.org_chrs[k['id_a']];
  chrs_b = PR.org_chrs[k['id_b']];

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

  #T = dict([T.items()[0]]);

  tests = 0;
            # Chromosome 1, chromosome2, index, score, pvalue)
  #start_clusts = [ (k[0], k[1], -1, None, None) for k in T.keys() ];
  start_clusts = [ (1, 5, 8221), (1, 5, 8480), (12, 2, 1521) ];
  sig_clusts   = [];
  C            = [];

  nd = cluster.null_dist(dist=null.cluster_null_score_strict_smart,
                 tests = 1,
                 scores=[ h[12] for h in hits],
                 hits=hits, hitA=hitA, hitB=hitB,
                 sperm=100000,
                 alpha = 0.05,
                 storage='null_dist_store_cmp_clt.dat');

  istack = [ (c1, c2, -1) for (c1,c2) in T.keys() ];
  #istack = [ (1, 5, 8221), (1, 5, 8480), (12, 2, 1521) ];
  #istack = [(1, 5, 8480)];
  istack = [ (4, 1, -1)];
  values = [];

  while len(istack) > 0:
    (c1, c2, index) = istack.pop();
    if index == None:
      continue;
    #fi

    t               = T[(c1,c2)];
    I, J, dist, ids = t[index];

    score, ex_a, ex_b, ue_a, ue_b, nz_ue_a, nz_ue_b = cluster.score_dendrogram_node(t, index, hits, chrs_a, chrs_b, hitA, hitB);

    p_perm = nd.pvalue_perm(score, len(ids), len(nz_ue_a), len(nz_ue_b));
    p_clt  = nd.pvalue_clt(score, len(ids), len(nz_ue_a), len(nz_ue_b));
    istack = [(c1, c2, I), (c1, c2, J)] + istack;
    values.append([(c1,c2,index),(len(ids), len(nz_ue_a), len(nz_ue_b)),p_perm, p_clt]);
    #values.append(p_clt);

    #if len(values) > 400:
    #  break;
    #fi

  #ewhile

  V = [ v for v in values if (v[3][1] != 0.0) and (v[2][1] != 0.0) ]

  Y     = np.array([ v[2][1] for v in V ], dtype=float);
  F     = np.array([ v[3][1] for v in V ], dtype=float);

  Ylog  = np.array([ math.log(v) for v in Y], dtype=float);
  Flog  = np.array([ math.log(v) for v in F], dtype=float);
  R2_v  = R2(Y,    F);
  R2log = R2(Ylog, Flog);
  rp    = np.corrcoef(F, Y);
  rplog = np.corrcoef(Flog, Ylog);

  return (V, R2_v, rp, R2log, rplog);
#edef

###############################################################################

def plot_clt_v_perm(V, thresh=0):

  V = [ v for v in V if (v[3][1] != 0.0) and (v[2][1] != 0.0) ];
  V = [ v for v in V if all(np.array(v[1], dtype=int) >= thresh) ];

  Y     = np.array([ v[2][1] for v in V ], dtype=float); # Perms
  F     = np.array([ v[3][1] for v in V ], dtype=float); # CLTs
  S     = np.array([ v[1][0] for v in V ], dtype=int);   # SIZE

  Ylog  = np.array([ -math.log(v) for v in Y], dtype=float);
  Flog  = np.array([ -math.log(v) for v in F], dtype=float);

  rp    = np.corrcoef(F, Y);
  rplog = np.corrcoef(Flog, Ylog);

  plt.cla();
  plt.plot([0,1], [0,1], color='r', linestyle='-', linewidth=2, zorder=1);
  plt.scatter(Y, F, zorder=2);
  plt.xlabel('Permutation p-values');
  plt.ylabel('CLT p-values');
  plt.axis([0,1,0,1]);
  plt.text(0.2,0.6,'r = %f' % rp[0][1]);
  plt.savefig('perms_v_clts_t%d.png' % thresh);
  plt.savefig('perms_v_clts_t%d.eps' % thresh)

  maxYlog = math.ceil(max(Ylog));
  maxFlog = math.ceil(max(Flog));
  plt.cla();
  plt.plot([0,max(maxYlog,maxFlog)], [0,max(maxYlog,maxFlog)], color='r', linestyle='-', linewidth=2, zorder=1);
  plt.scatter(Ylog, Flog, zorder=2);
  plt.xlabel('Permutation p-values (-log)');
  plt.ylabel('CLT p-values (-log)');
  plt.axis([0,maxYlog,0,maxFlog]);
  plt.text(2, 20, 'r = %f' % rplog[0][1]);
  plt.savefig('perms_v_clts_t%d_log.png' % thresh);
  plt.savefig('perms_v_clts_t%d_log.eps' % thresh);

  return (thresh, len(V), rp[0][1], rplog[0][1]);
#edef

###############################################################################

def R2(Y, F):
  muY = np.mean(Y);
  ssres = sum([ (y - f)**2 for (y,f) in zip(Y,F) ]);
  sstot = sum([ (y - muY)**2 for y in Y ]);
  R2 = 1 - ssres/sstot;
  return R2;
#edef

###############################################################################

def peaks_analysis_reconstruct_tree(PR, k):

  Cs = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd'])];
  C = {};
  T = {};

  # Split by chrs
  for c in Cs:
    key = (c[0],c[3]);
    if not(key in C.keys()):
      C[key] = [];
    #fi
    C[key].append(c);
  #efor

  for key in C.keys():
    clusters = C[key];
    t = [];
    S = [];

    clusters = sorted(clusters, key=(lambda x: x[11]), reverse=False);
    for (i,clust) in enumerate(clusters):
      ovs      = [ (clust[1] <= x[1] and clust[2] >= x[2]) for (ix,x) in S ];
      children = [ ix for (o,(ix,s)) in zip(ovs,S) if o ];
      S        = [ s for (o,s) in zip(ovs,S) if not(o) ] + [ (i,clust) ];
      t.append((clust[7][0],tuple(children)));
    #efor
    root_nodes = [ i for (i,c) in S ];
    T[key] = (clusters, t, root_nodes);
  #efor

  return T;
#efor

###############################################################################

def peaks_prominance(T, root=-1, minimum=True):

  if minimum:
    T = [ (-h,c) for (h,c) in T ];
  #fi

  H,Children = zip(*T);

  istack     = [ root ];
  Hstack     = [ ];
  Visited    = [ False for i in xrange(len(T)) ];
  prominence = [ None  for i in xrange(len(T)) ];
  right_cols = [ None  for i in xrange(len(T)) ];

  while len(istack) > 0:
    index  = istack[-1];
    height = H[index];
    childs = Children[index];
    visit  = Visited[index];

    if visit:
      for i in xrange(len(Hstack)-1, -1, -1):
        h = Hstack[i];
        if h > height:
          break;
        #fi
      #efor
      leftcol  = min(Hstack[i:]);

      if len(childs) == 0:
        rightcol      = [min(0,height)];
        rightcol_path = rightcol;
      elif max([H[c] for c in childs]) > height:
        rightcol      = [height];
        rightcol_path = rightcol;
      else:
        rightpaths     = [ right_cols[c] for c in childs ];
        node_rightcols = [];
        for rp in rightpaths:
          for (i,h) in enumerate(rp):
            if h > height:
              break;
            #fi
          #efor
          node_rightcols.append(min(rp[0:i] + [0]));
        #efor
        rightcol_path = [ height ] +  rightpaths[node_rightcols.index(max(node_rightcols))];
        rightcol = node_rightcols;
      #fi

      prominence[index] = height - max(leftcol, max(rightcol));
      right_cols[index] = rightcol_path;

      istack.pop();
      Hstack.pop();
    else:
      Visited[index] = True;
      Hstack.append(height);
      for child in childs:
        istack.append(child);
      #efor
    #fi

  #ewhile

  return prominence;

#edef
    
###############################################################################

def peaks_analysis(PR, k, N=15):

  NK = {};
  for key in k.keys():
    NK[key] = k[key];
  #efor
  NK['cut'] = 'peaks_%d' % N;

  trees = peaks_analysis_reconstruct_tree(PR, k);
  PROM  = {};
  SCP   = [];

  for key in trees.keys():
    C, T, root_nodes = trees[key];
    P = [ None for i in xrange(len(T)) ];
    for rn in root_nodes:
      prom = peaks_prominance(T, rn, minimum=True);
      for i in xrange(len(T)):
        if prom[i] != None:
          P[i] = prom[i];
        #fi
      #efor
    #efor
    PROM[key] = P;
    SCP       = SCP + [ C[i] for (i,s) in sorted(enumerate(P), key=lambda x:x[1])[-N:] ];
      
    #efor
  #efor
  PR.hit_clusters[(NK['id_a'], NK['id_b'], NK['linkage_type'], NK['alpha'], NK['cut'], NK['nd'])] = SCP; 
  
  return NK;

#edef

###############################################################################

def cmp_runs(PRs, ks, names, prefix):

  feats = [];

  for (PR, k, name) in zip(PRs, ks, names):
    fs = [];
    clusts = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd'])];

    for c in clusts:
      feat_a_size   = c[2] - c[1];
      feat_b_size   = c[5] - c[4];
      feat_rel_size = float(feat_a_size) / float(feat_b_size);
      feat_max_size = max(feat_a_size, feat_b_size);

      feat_n     = c[7][1];
      feat_nue_a = c[7][2];
      feat_nue_b = c[7][3];
      feat_nue   = feat_nue_a + feat_nue_b;
      feat_ratio = float(feat_n+1) / float(feat_nue_a + feat_nue_b+1);

      feat_score  = c[7][0];
      feat_pvalue = c[8];

      feat_ngenes_a    = len(c[9]);
      feat_ngenes_b    = len(c[10]);
      feat_gene_den    = float(feat_ngenes_a+1) / float(feat_ngenes_b+1);
      feat_total_genes = feat_ngenes_a + feat_ngenes_b;

      f = (feat_a_size, feat_b_size, feat_rel_size, feat_max_size,
           feat_n, feat_nue_a, feat_nue_b, feat_nue, feat_ratio,
           feat_score, feat_pvalue,
           feat_ngenes_a, feat_ngenes_b, feat_gene_den, feat_total_genes);
      fs.append(f);
    #efor
    feats.append(zip(*fs));
  #efor

  plot_feats(feats, 8, 9,
             "Cluster conservation ratio vs. cluster score",
             "Cluster conservation ratio",
             "Cluster score",
             names,
             "%s_cluster_ratio_vs_score" % prefix, logx=True);

  plot_feats(feats, 9, 10,
             "Cluster score vs. Cluster p-value",
             "Cluster score",
             "Cluster p-value",
             names,
             "%s_cluster_score_vs_pvalue" % prefix, logy=True, YMIN=10**-20, YMAX=1);

  plot_feats(feats, 11, 12,
             "Number of genes in each organism",
             "# of genes in org1",
             "# of genes in org2",
             names,
             "%s_n_genes" % prefix, log=True);

  plot_feats(feats, 9, 13,
             "Cluster score vs. ratio of # genes in org1/org2",
             "Cluster score",
             "# of genes in org1/org2",
             names,
             "%s_cluster_score_vs_density" % prefix, logy=True);

  plot_feats(feats, 0, 1,
             "Cluster sizes on each organism",
             "Cluster size on org 1",
             "Cluster size on org 2",
             names,
             "%s_cluster_sizes" % prefix, legendloc="lower right", log=True);

  plot_feats(feats, 3, 4,
             "Cluster size vs. number of hits in cluster",
             "Cluster size",
             "Number of hits in cluster",
             names,
             "%s_clust_size_vs_hits" % prefix, log=True, legendloc="upper left");

  plot_feats(feats, 3, 14,
             "Cluster size vs. number of genes in cluster",
             "Cluster size",
             "Number of genes in cluster",
             names,
             "%s_clust_size_vs_genes" % prefix, log=True, legendloc="upper left")

  plot_feats(feats, 3, 8,
             "Cluster size vs conservation_ratio",
             "Cluster size",
             "Conservation ratio",
             names,
             "%s_clust_size_vs_ratio" % prefix, log=True);

  plot_feats(feats, 4, 7,
             "Number of hits in cluster vs unaccounted exons",
             "Number of hits",
             "Total number of unaccounted exons",
             names,
             "%s_hits_vs_nue_total" % prefix, log=True);

  return feats;
#edef 

###############################################################################

def plot_feats(feats, i1, i2, title, xlab, ylab, labels, out, log=False, logx=False, logy=False, legendloc=4, XMIN=None, XMAX=None, YMIN=None, YMAX=None):
  colors = [ 'b', 'g', 'r', 'c', 'm', 'y', 'k' ];
  plt.cla();
  xmin  = float("inf");
  xmax  = float("-inf");
  ymin  = float("inf");
  ymax  = float("-inf");
  for (F,c, l) in zip(feats, colors[0:len(feats)], labels):
    plt.scatter(F[i1], F[i2], c=c, label=l);
    xmin  = min(xmin,  min(F[i1]));
    xmax  = max(xmax,  max(F[i1]));
    ymin  = min(ymin,  min(F[i2]));
    ymax  = max(ymax,  max(F[i2]));
  #efor
  plt.xlabel(xlab);
  plt.ylabel(ylab);
  plt.title(title);
  if logx or log:
    plt.xscale('log');
    if xmin == 0:
      xmin = 1;
    #fi
  #fi
  if logy or log:
    plt.yscale('log');
    if ymin == 0:
      ymin = 1;
    #fi
  #fi

  xmin = xmin if (XMIN == None) else XMIN;
  xmax = xmax if (XMAX == None) else XMAX;
  ymin = ymin if (YMIN == None) else YMIN;
  ymax = ymax if (YMAX == None) else YMAX;

  plt.xlim([xmin, xmax]);
  plt.ylim([ymin, ymax]);

  print (xmin, xmax, ymin, ymax);

  plt.legend(loc=legendloc);

  plt.savefig('%s.png' % out);
  plt.savefig('%s.eps' % out);
#edef

###############################################################################

def cmp_run_overlaps(PRa, ka, PRb, kb):

  clusts_a = PRa.hit_clusters[(ka['id_a'], ka['id_b'], ka['linkage_type'], ka['alpha'], ka['cut'], ka['nd'])];
  clusts_b = PRb.hit_clusters[(kb['id_a'], kb['id_b'], kb['linkage_type'], kb['alpha'], kb['cut'], kb['nd'])];

  chr_groups_a = {};

  total_coverage_a = 0;
  total_coverage_b = 0;
  overlap          = 0;

  for c in clusts_a:
    k = (str(c[0]), str(c[3]));
    if k not in chr_groups_a:
      chr_groups_a[k] = [];
    #fi
    chr_groups_a[k].append(c);
    total_coverage_a = total_coverage_a + (c[2] - c[1]) + (c[5] - c[4]);
  #efor

  for c in clusts_b:
    k = (str(c[0]), str(c[3]));
    total_coverage_b = total_coverage_b + (c[2] - c[1]) + (c[5] - c[4]);

    if k not in chr_groups_a:
      continue;
    #fi

    rel_clusts = chr_groups_a[k];
    for rc in rel_clusts:
      ov_a = util.overlap((c[1], c[2]), (rc[1], rc[2]));
      ov_b = util.overlap((c[4], c[5]), (rc[4], rc[5]));
      if ov_a > 0 and ov_b > 0:
        print k
        print (c[1], c[2]), (rc[1], rc[2]), "->", ov_a
        print (c[4], c[5]), (rc[4], rc[5]), "->", ov_b
        overlap = overlap + ov_a + ov_b;
      #fi
    #efor

  print overlap, total_coverage_a, total_coverage_b

  return (float(overlap) / float(total_coverage_a + total_coverage_b));

#edef

###############################################################################
