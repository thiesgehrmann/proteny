#!/usr/bin/python

import sys;
sys.path.append('utils');

import circos_data;
import proteny as pp;
import data;

###############################################################################

savename = 'PR.proteny';
circdir  = 'circos';

###############################################################################
# DATA

  # Download data for Schizophyllum commune
sc = data.schco2();
  # Download data for Agaricus bisporus
ag = data.agabi();

  # Prepare the protenty analysis
PR = pp.proteny();

  # Add our organisms.
id_a = PR.add_org(*sc, isfile=False);
id_b = PR.add_org(*ag, isfile=False);

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
cdata = circos_data.circos_data();
cdata.write_data(PR, '%s/data' % circdir);


