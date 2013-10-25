import sys;
sys.path.append('./utils');

import proteny as ps;
import util;

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

def hit_cluster_overlap(PR, k):
  """ Identify clusters which are overlapped"""

  if ((k['id_a'], k['id_b'], k['linkage_type'], k['alpha']) not in PR.hit_clusters):
    print "You must run hit_cluster() first!";
    return None;
  #fi

  C = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'])];

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


  C  = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'])];
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
  C  = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'])][ci];

  #('bri1',   [ ('schco2', 'scaffold_1',  5346259, 5348750), ('agabi', 'scaffold_1',  1865954, 1869064) ] ),

  reg = (name if not(name == None) else str(ci),
         [ ( PR.org_names[k['id_a']], C[0], C[1], C[2] ),
           ( PR.org_names[k['id_b']], C[3], C[4], C[5] ) ] );

  return reg;
#edef

###############################################################################
