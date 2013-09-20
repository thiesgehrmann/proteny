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
human = data.human();
  # Download data Aniger 513_88
mouse = data.mouse();

  # Prepare the protenty analysis
PR = pp.proteny();

  # Add our organisms.
id_a = PR.add_org(*an1, isfile=False);
id_b = PR.add_org(*an2, isfile=False);

###############################################################################
# ANALYSIS

  # Run BLAST
k = PR.blast(id_a=id_a, id_b=id_b);
  # Calculate window scores for these results
k = PR.windows(k);
  # Calculate genomic distances between ALL hits on different chromosomes
k = PR.hit_distance(k);
  # Using a specified linkage, perform agglomerative clustering
k = PR.cluster_linkage(k);
  # Cut the dendrogram.
k = PR.cluster_hits(k);

###############################################################################
# SAVE

  # Save our calculations
PR.save(savename);

###############################################################################
# VISUALIZE

  # Export data into a CIRCOS format
files = cdata.write_data(PR, '%s/data' % circdir);

  # Visualize relationships between chromosomes
anchrs = [ "human_%d" % (i+1)  for i in xrange(an2[2].Shape()()) ];
for i in xrange(an1[2].Shape()()):
 cc.circos_chr(files, anchrs + ["mouse_%d" % (i+1)], ["mouse%d=0.4r" % (i+1)], circdir, "mouse_%d_view" % (i+1) );
#efor
