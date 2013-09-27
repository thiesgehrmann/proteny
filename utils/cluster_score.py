# Need:
# cluster chr

import bisect;

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
