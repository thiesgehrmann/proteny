from ibidas import *;

###############################################################################

def reindex_blast(br):
  nbr = br.Shape()();
  ind = [i for i in xrange(nbr)];

  F = br.Get(_.a_chrid, _.a_geneid, _.a_exonid, \
             _.b_chrid, _.b_geneid, _.b_exonid, \
             _.a_start,                         \
             _.a_start,                         \
             _.b_start,                         \
             _.b_start,                         \
             _.pident,                          \
             _.evalue,                          \
             _.bitscore);
  F = Rep(zip(ind, *F()))
  F = F / ('i', 'a_chrid', 'a_geneid', 'a_exonid', 'b_chrid', 'b_geneid', 'b_exonid', 'a_start', 'a_end', 'b_start', 'b_end', 'pident', 'evalue', 'bitscore');

  return F.Copy();
#edef

###############################################################################

def sort_org(br, hits):
  chrs = br.chrid.Unique().Sort()();

  chr_hits = [];

  for (i, chr) in enumerate(chrs):
    chr_br = br[_.chrid == chr].Sort(_.start);
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

  overlap = reduce(lambda x,y: x ^ y, ids);

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
###############################################################################
###############################################################################
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

  a_chr   = hits[0][1];
  a_start = min([ h[2] for h in hits_a]);
  a_end   = max([h[3] for h in hits_a]);
  b_chr   = hits[0][4];
  b_start = min([ h[2] for h in hits_b]);
  b_end   = max([h[3] for h in hits_b]);
  n_hits  = len(hits);
  score   = sum(hits_scores);

  return (a_chr, hits[0][0][0][0], a_start, a_end, b_chr, hits[0][0][1][0], b_start, b_end, n_hits, score, prots_a, prots_b);
#edef;


