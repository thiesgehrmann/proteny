#!/usr/bin/python

import sys;
sys.path.append('utils');
sys.path.append('visualization');

import circos_data as cdata;
import circos_region as cr;
import circos_chr    as cc;
import proteny as ps;
import data;
import cluster_null as null;
import postprocessing as pp;

###############################################################################

savename = 'PR_aniger.proteny';
circdir  = 'visualizations/aniger';

###############################################################################
# DATA

  # Download data for Aniger n402
an1 = data.aniger_n402();
  # Download data Aniger 513_88
an2 = data.aniger_513_88();

  # Prepare the protenty analysis
PR = ps.proteny();

  # Add our organisms.
id_a = PR.add_org(*an1, isfile=False);
id_b = PR.add_org(*an2, isfile=False);

###############################################################################
# ANALYSIS

#k = PR.analyze(id_a=id_a, id_b=id_b);
k  = PR.blast(id_a=id_a, id_b=id_b);
k  = PR.hit_index(k);
k  = PR.hit_distance(k);
k  = PR.hit_dendrogram(k);
k1 = PR.hit_cluster(k, cut='simple', nd=null.cluster_null_score_strict);
k2 = PR.hit_cluster(k, cut='greedy', nd=null.cluster_null_score_strict);


###############################################################################
# SAVE

  # Save our calculations
PR.save(savename);

###############################################################################
# VISUALIZE

def viz(PR, k, circdir):
  files = cdata.write_data(PR, k, '%s/data' % circdir);

  KS = PR.key_s(k);

  regions = [ # ('strange',  [ ('513.88', '1',  0, 100), ('513.88', '8', 0, 100), ('n402', '4', 0, 100) ] ), 
              ('strange2', [ ('513.88', '12', 0, 100), ('513.88', '5', 0, 100), ('n402', '5', 0, 100) ] ) ];

  anchrs = [ "513.88_%d" % (i+1)  for i in xrange(an2[2].Shape()()) ];
  for i in xrange(an1[2].Shape()()):
   cc.circos_chr(files, anchrs + ["n402_%d" % (i+1)], ["n402_%d=0.4r" % (i+1)], circdir, "scaffold_%02d_%s" % (i+1, KS) );
  #efor

    # Visualize a region
  for reg in regions:
    cr.circos_region(files, reg, 30000, circdir, ('RegKnown_' + KS + '_' + reg[0]));
    print 'RegKnown_' + KS + '_' + reg[0];
  #efor
 
#edef


  # Produce circos visualizations
viz(PR, k1, 'visualizations/basid_more/greedy');
viz(PR, k2, 'visualizations/basid_more/simple');

###############################################################################
###############################################################################

#IR  = pp.hit_cluster_overlap(PR, k);
#IRs = [ ir for ir in IR if ir[0] == '18' ];

#IRs = [ ('correct_mapping', [ ('n402', '21', 884643, 917655), ('513.88', '8', 880417, 909428), ('513.88', '14', 482348, 497064) ]) ];
#for reg in IRs:
#  cr.circos_region(files, reg, 100, circdir, ('RegOverlapping_' + KS + '_' + reg[0]));
#  print 'RegOverlapping_' + KS + '_' + reg[0];
##efor


  # Visualize relationships between chromosomes
anchrs = [ "513.88_%d" % (i+1)  for i in xrange(an2[2].Shape()()) ];
for i in xrange(an1[2].Shape()()):
 cc.circos_chr(files, anchrs + ["n402_%d" % (i+1)], ["n402_%d=0.4r" % (i+1)], circdir, "scaffold_%02d_%s" % (i+1, KS) );
#efor

