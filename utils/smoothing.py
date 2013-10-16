import numpy as np;

from ibidas.utils.util import debug_here;
from ibidas import *;

###############################################################################

def sort_org(br, hits):
  """ O, H = sort_org(br, hits)
      br:   The output of reindex_blast
      hits: A list which contains a list of information about a hit, prefixed by an empty list.

      Output:
        O: A list of hits sorted by start location
        H: The updated list hits, containing indexes to find the sorted element
  """

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
  """L = get_reg_hits(H, O, k, ws)
     H:  The output list of sort_org
     O:  The output list of sort_org
     k:  The ID of the hit
     ws: The size of the window around the hit.

     Outputs:
      L: A list of hit IDs present in the window.
  """

  hit = H[k];
  ids = [];

  for (i, o) in enumerate(O):
    reg = hit[0][i];
    ids.append(set(get_reg_window(o[reg[0]], reg[1], ws)));
  #efor

  ovl = reduce(lambda x,y: x & y, ids);

  ret = [];
  for i in ovl:
    ret.append(H[i]);
  #efor

  return ret;
#def

###############################################################################

def get_reg_window(O, k, ws):
  """ids = get_reg_window(O, k, ws):
     O:  The output list of sort_org
     k:  The ID if the hit
     ws: The size of the window around the hit

     Outputs:
       ids: A list of IDs of hits around the hit.
  """

  rhit = O[k];
  olen = len(O);

  min_start = rhit[2] - ws;
  max_end   = rhit[3] + ws;

  dstream_i = k;
  ustream_i = k;

  #debug_here();

  while dstream_i-1 > 0:
    if O[dstream_i-1][2] > min_start:
      dstream_i -= 1;
    else:
      break;
    #fi
  #ewhile
  while ustream_i+1 < olen:
    if O[ustream_i+1][3] < max_end:
      ustream_i += 1;
    else:
      break;
    #fi
  #ewhile

  return [ o[0] for o in O[dstream_i:ustream_i+1] ];
#edef

###############################################################################

