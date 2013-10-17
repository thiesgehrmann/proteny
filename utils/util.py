from ibidas import *;

import numpy as np;

###############################################################################

def chr_pair_group(hits):
  """HC = chr_pair_group(H)
     H: Output list from sort_org

     Group hits by chromosomes they exist on
     NOTE: DOES NOT PRESERVE INDICES OF O LISTS!
     Outputs:
       HC: A dictionary of hits per pair of chromosomes
  """

  H_chrs = {};

  for (i, hit) in enumerate(hits):
    a_chr1 = hit[1];
    a_chr2 = hit[4];
    k = (a_chr1, a_chr2);
    if k not in H_chrs:
      H_chrs[k] = [];
    #fi
    H_chrs[k].append(i);
  #efor

  return H_chrs;

#edef

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

def reindex_blast(br):
  """ B + reindex_blast(br)
      br: A rep of BLAST results.
          Must contain at least:
          _.a_chrid, _.a_geneid, _.a_exonid
          _.b_chrid, _.b_geneid, _.b_exonid
          _.a_start, _.a_end
          _.b_start, _.b_end
          _.pident, _.evalue, _.bitscore

      Output:
        B: A Rep of blast results with an index
  """

  nbr = br.Shape()();
  ind = [i for i in xrange(nbr)];

  F = br.Get(_.a_chrid, _.a_geneid, _.a_exonid, \
             _.b_chrid, _.b_geneid, _.b_exonid, \
             _.a_start,                         \
             _.a_end,                           \
             _.b_start,                         \
             _.b_end,                           \
             _.pident,                          \
             _.evalue,                          \
             _.bitscore);
  F = Rep(zip(ind, *F()))
  F = F / ('i', 'a_chrid', 'a_geneid', 'a_exonid', 'b_chrid', 'b_geneid', 'b_exonid', 'a_start', 'a_end', 'b_start', 'b_end', 'pident', 'evalue', 'bitscore');

  return F.Copy();
#edef

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

