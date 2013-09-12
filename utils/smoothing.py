from ibidas import *;

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
