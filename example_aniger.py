#!/usr/bin/python

import sys;

from visualization import circos_data as cdata;
from visualization import circos_region as cr;
import circos_chr    as cc;
import proteny as ps;
import data;
from utils import cluster_null as null;
import postprocessing as pp;

###############################################################################

savename = 'PROTENY_OUTPUT_aspni.proteny';
circdir  = 'PROTENY_OUTPUT/aspni';


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

  # Run analysis
k = PR.analyze(id_a=id_a, id_b=id_b, cut='deeper_greater', nd=null.cluster_null_score_strict_smart, alpha=0.05, ngenes_threshold=2, conservation_ratio=2);
  # Save! It's important!
PR.save(savename);

  # Produce circos visualizations
viz(PR, k, circdir);



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
viz(PR, k, 'visualizations/basid_more/greedy');

###############################################################################
###############################################################################

