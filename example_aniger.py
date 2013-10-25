#!/usr/bin/python

import sys;
sys.path.append('utils');
sys.path.append('visualization');

import circos_data as cdata;
import circos_region as cr;
import circos_chr    as cc;
import proteny as pp;
import data;

###############################################################################

savename = 'PR_aniger.proteny';
circdir  = 'circos';

###############################################################################
# DATA

  # Download data for Aniger n402
an1 = data.aniger_n402();
  # Download data Aniger 513_88
an2 = data.aniger_513_88();

  # Prepare the protenty analysis
PR = pp.proteny();

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

           #('bri1',   [ ('schco2', 'scaffold_1',    5346259, 5348750), ('agabi', 'scaffold_1',  1865954, 1869064) ] ),
regions = [ ('strange',[ ('aniger_n402',   'ACJE01000001',  436263,  466251),
                         ('aniger_513_88', 'An01',         1618193, 1644408),
                         ('aniger_513_88', 'An11',         1443277, 1470326) ]) ];

  # Visualize a region
for reg in regions:
  cr.circos_region(files, reg, 30000, circdir, ('aniger_reg_%s' % reg[0]));
#efor


  # Visualize relationships between chromosomes
anchrs = [ "aniger_513_88_An%02d" % (i+1)  for i in xrange(an2[2].Shape()()) ];
for i in xrange(an1[2].Shape()()):
 cc.circos_chr(files, anchrs + ["aniger_n402_ACJE010000%02d" % (i+1)], ["aniger_n402_ACJE010000%02d=0.4r" % (i+1)], circdir, "An%02d" % (i+1) );
#efor

