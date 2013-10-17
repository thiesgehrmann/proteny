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

savename = 'PR_human12_mouse6.proteny';
circdir  = 'circos';

###############################################################################
# DATA

  # Download data for Aniger n402
#human = data.human();
human = data.example(data.human, '12')
  # Download data Aniger 513_88
#mouse = data.mouse();
mouse = data.example(data.mouse, '6')

  # Prepare the protenty analysis
PR = pp.proteny();

  # Add our organisms.
id_a = PR.add_org(*human, isfile=False);
id_b = PR.add_org(*mouse, isfile=False);

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
files = cdata.write_data(PR, '%s/data' % circdir);

  # Visualize relationships between chromosomes
anchrs = [ "human_%d" % (i+1)  for i in xrange(human[2].Shape()()) ];
for i in xrange(mouse[2].Shape()()):
 cc.circos_chr(files, anchrs + ["mouse_%d" % (i+1)], ["mouse%d=0.4r" % (i+1)], circdir, "mouse_%d_view" % (i+1) );
#efor

