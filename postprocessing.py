import sys;
sys.path.append('./utils');

import proteny as ps;
import util;

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
