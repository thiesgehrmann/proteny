
###############################################################################

def reindex_blast(br):
  nbr = br.Shape()();
  ind = [i for i in xrange(nbr)];

  F = br.Get(_.a_assemblyid, _.a_proteinid, _.a_exonid, \
             _.b_assemblyid, _.b_proteinid, _.b_exonid, \
             (_.a_start + (_.qstart * 3 )),                         \
             (_.a_start + (_.qend * 3 )),                           \
             (_.b_start + (_.sstart * 3 )),                         \
             (_.b_start + (_.send * 3 )),                           \
             _.pident,                                              \
             _.evalue,                                              \
             _.bitscore);
  F = Rep(zip(ind, *F()))
  F = F / ('i', 'a_assemblyid', 'a_proteinid', 'a_exonid', 'b_assemblyid', 'b_proteinid', 'b_exonid', 'a_start', 'a_end', 'b_start', 'b_end', 'pident', 'evalue', 'bitscore'); 

  return F.Copy();
#edef

###############################################################################

def sort_org(br, hits):
  chrs = br.assemblyid.Unique().Sort()();

  chr_hits = [];

  for (i, chr) in enumerate(chrs):
    chr_br = br[_.assemblyid == chr].Sort(_.start);
    chr_hits.append([]);
    for (j, hit) in enumerate(zip(*chr_br())):
      hit_id, hit_chr, hit_start, hit_end = hit;
      chr_hits[i].append(hit);
      hits[hit_id][0].append((i,j));
    #efor
  #efor
  return (chr_hits, hits);
#edef

###############################################################################

def get_reg_hits(H, O, k, ws):
  hit = H[k];
  ids = [];

  for (i, o) in enumerate(O):
    reg = hit[0][i];
    ids.append(set(get_reg_window(o[reg[0]], reg[1], ws)));
  #efor

  overlap = reduce(lambda x,y, x ^ y, ids);

  ret = [];
  for i in overlap:
    ret.append(H[i]);
  #efor

  return ret;
#def

###############################################################################

def get_reg_window(O, k, ws):
  rhit = O[k];
  olen = len(O);

  min_start = rhit[2] - ws;
  max_end   = rhit[3] + ws;

  dstream_i = k;
  ustream_i = k;

  while dstream_i > 0:
    if O[dstream_i][2] > min_start:
      dstream_i -= 1;
    else:
      break;
    #fi
  #ewhile
  while ustream_i < olen:
    if O[ustream_i][3] < max_end:
      ustream_i += 1;
    else:
      break;
    #fi
  #ewhile

  return [ o[0] for o in O[dstream_i:ustream_i+1] ];
#edef

###############################################################################

def score(L):
  total = sum([l[-1] for l in L]) / len(L);
  return total;
#efor

###############################################################################

def overlap(r1, r2):
  """
    d = overlap(r1, r2)

      Returns the overlap d of two regions specified.
      If the regions do not overlap, then the negative distance between them is returned.

      r1 = (start1, end1)
      r2 = (start2, end2);

    s1      e1
    |       |
    ####1####
           ||  <---------  (overlap)
           ####2####
           |       |
           s2      e2
  """

  if r1[0] < r2[0]:
    e1 = r1[1];
    s2 = r2[0];
  else:
    e1 = r2[1];
    s2 = r1[0];
  #fi

  ov = e1 - s2;

  return ov;
#edef

###############################################################################

def distance(H, O, i, j):
  h1 = H[i];
  h2 = H[j];

    # Different chromosomes
  if (h1[1] != h2[1]) or (h1[4] != h2[4]):
    return float("inf");
  #fi

    # Same exons
  if (h1[2] == h2[2]) and (h1[5] == h2[5]) and (h1[3] == h2[3]) and (h1[6] == h2[6]):
    return 0;
  #fi

    # Oh Jesus save me.
  h1_a = O[0][h1[0][0][0]][h1[0][0][1]];
  h1_b = O[1][h1[0][1][0]][h1[0][1][1]];
  h2_a = O[0][h2[0][0][0]][h2[0][0][1]];
  h2_b = O[1][h2[0][1][0]][h2[0][1][1]];
  ov1 = overlap(h1_a[-2:], h2_a[-2:]);
  ov2 = overlap(h1_b[-2:], h2_b[-2:]);

    # If they all overlap, distance is 0.
  if ov1 > 0 and ov2 > 0:
    return 0;
  #fi

  return -(min(0, ov1) + min(0, ov2));
#edef

###############################################################################

def condenseddm(H, O, L):
  nl = len(L);

  cdm = [];

  for i in xrange(nl-1):
    for j in xrange(i+1, nl):
      cdm.append(distance(H, O, L[i], L[j]));
    #efor
  #efor

  return cdm;
#edef

###############################################################################

def chr_pair_group(H):
  """ NOTE: DOES NOT PRESERVE INDICES OF O LISTS! """

  H_chrs = {};

  for i in xrange(len(H)):
    a_chr1 = H[i][1];
    a_chr2 = H[i][4];
    k = (a_chr1, a_chr2);
    if k not in H_chrs:
      H_chrs[k] = [];
    #fi
    H_chrs[k].append(i);
  #efor

  return H_chrs;
#edef

###############################################################################

def clust_description(H, O, scores, C):
  hits        = [ H[i] for i in C ];
  hits_scores = [ scores[i][-4] for i in C ];
  hits_a      = [ O[0][h[0][0][0]][h[0][0][1]] for h in hits ];
  hits_b      = [ O[1][h[0][1][0]][h[0][1][1]] for h in hits ];
  prots_a     = set([ h[2] for h in hits ]);
  prots_b     = set([ h[5] for h in hits ]);

  a_scaff = hits[0][1];
  a_start = min([ h[2] for h in hits_a]);
  a_end   = max([h[3] for h in hits_a]);
  b_scaff = hits[0][4];
  b_start = min([ h[2] for h in hits_b]);
  b_end   = max([h[3] for h in hits_b]);
  n_hits  = len(hits);
  score   = sum(hits_scores);

  return (a_scaff, hits[0][0][0][0], a_start, a_end, b_scaff, hits[0][0][1][0], b_start, b_end, n_hits, score, prots_a, prots_b);
#edef;

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

BR = Load("results.Schco.Agabi.blast");
F  = reindex_blast(BR);

hits = F.Get(_.i, _.a_assemblyid, _.a_proteinid, _.a_exonid, _.b_assemblyid, _.b_proteinid, _.b_exonid, _.pident, _.evalue, _.bitscore);
BR_a = F.Get(_.i, _.a_assemblyid, _.a_start, _.a_end) / ('i', 'assemblyid', 'start', 'end');
BR_b = F.Get(_.i, _.b_assemblyid, _.b_start, _.b_end) / ('i', 'assemblyid', 'start', 'end');

brs  = [ BR_a, BR_b ];

H = [ [[]] + list(x[1:]) for x in zip(*hits()) ]
O = [];

for br in brs:
  org_ind, H = sort_org(br, H);
  O.append(org_ind);
#efor

RS = [];
ws = [ 200, 2000, 20000, 40000, 100000 ];
hlen = len(H);

for w in ws:
  R = [];
  for h in xrange(hlen):
    print "%d, %d / %d" % (w, h, hlen);
    R.append(get_reg_hits(H, O, h, w));
  #efor
  RS.append(R);
#efor

scores = [[ score(RS[i][k]) for i in xrange(5) ] for k in xrange(hlen)]

k = sorted(enumerate(scores), key=lambda x: x[1][2]);

###############################################################################

from scipy.cluster import hierarchy;
import numpy as np;

HC = H_chrs_group(H);
metric = lambda i, j: distance(H, O, int(i[0]), int(j[0]));

X      = np.array([ [x] for x in HC[(u'scaffold_10', u'scaffold_10')]], dtype=np.dtype('u8'));
hclust = hierarchy.fclusterdata(X, 0, metric=metric, criterion='distance');

#L = {};

#for k in HC.keys():
#  print k
#  hc      = HC[k];
#  if len(hc) < 2:
#    continue;
#  #fi
#
#  cdm     = np.array(condenseddm(H, O, hc), dtype=np.dtype('u8'));
#  if len(cdm) < 2:
#    continue;
#  #fi
#
#  L[k] = hierarchy.single(cdm);
##efor
#Save(L, 'linkage_single.dat');

clust_sizes = [ 200, 2000, 10000, 20000, 50000, 100000, 200000, 500000 ];
clust_sizes = [ 20000 ];

L = Load('linkage_single.dat');

for cs in clust_sizes:

  HC_clust = {};
  for k in L.keys():
    linkage = L[k];
    hc      = HC[k];
    hclust  = hierarchy.fcluster(linkage, cs, criterion='distance');
    nclust  = max(hclust);
    clusts  = [ [] for i in xrange(nclust) ];

    for (h, c) in zip(hc, hclust):
      clusts[c-1].append(h);
    #efor
    HC_clust[k] = clusts;
  #efor

  hit_clusts = [];
  for k in HC_clust.keys():
    clusts = HC_clust[k];
    for C in clusts:
      cd = clust_description(H, O, scores, C);
      hit_clusts.append(cd);
    #efor
  #efor

  CIRCOS_DIR = 'test_circos';
  names      = [ "Schco", "Agabi" ];

#  aprot = set([ 72717, 136031, 188164, 192455, 201136, 214600, 215236]);

  fd = open('%s/data/clusters_%d_avg_sum_subset.tsv' % (CIRCOS_DIR, cs), 'w');
  for (i, h) in enumerate(hit_clusts):
#    if len(aprot & h[-1]) == 0:
#      continue;
#    #fi
    fd.write('clusters_%s_%s_%010d\t%s_%s\t%d\t%d\ta_chrid=%d,nhits=%d,score=%f,prots_a=%s\n' % (names[0], names[1], i, names[0], h[0], h[2], h[3], h[1], h[8], h[9], ';'.join(['%s' % id for id in h[10]])));
    fd.write('clusters_%s_%s_%010d\t%s_%s\t%d\t%d\tb_chrid=%d,nhits=%d,score=%f,prots_b=%s\n' % (names[0], names[1], i, names[1], h[4], h[6], h[7], h[5], h[8], h[9], ';'.join(['%s' % id for id in h[11]])));
    i += 1;
  #efor
  fd.close();
#efor


###############################################################################



