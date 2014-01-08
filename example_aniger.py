#!/usr/bin/python

import sys;
sys.path.append('utils');
sys.path.append('visualization');

import circos_data as cdata;
import circos_region as cr;
import circos_chr    as cc;
import proteny as ps;
import data;

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

k = PR.analyze(id_a=id_a, id_b=id_b);

###############################################################################
# SAVE

  # Save our calculations
PR.save(savename);

###############################################################################
# VISUALIZE

  # Export data into a CIRCOS format
files = cdata.write_data(PR, k, '%s/data' % circdir);

KS = PR.key_s(k);

IR  = pp.hit_cluster_overlap(PR, k);
IRs = [ ir for ir in IR if ir[0] == '18' ];

IRs = [ ('correct_mapping', [ ('n402', '21', 884643, 917655), ('513.88', '8', 880417, 909428), ('513.88', '14', 482348, 497064) ]) ];
for reg in IRs:
  cr.circos_region(files, reg, 100, circdir, ('RegOverlapping_' + KS + '_' + reg[0]));
  print 'RegOverlapping_' + KS + '_' + reg[0];
#efor


  # Visualize relationships between chromosomes
anchrs = [ "513.88_%d" % (i+1)  for i in xrange(an2[2].Shape()()) ];
for i in xrange(an1[2].Shape()()):
 cc.circos_chr(files, anchrs + ["n402_%d" % (i+1)], ["n402_%d=0.4r" % (i+1)], circdir, "scaffold_%02d_%s" % (i+1, KS) );
#efor

