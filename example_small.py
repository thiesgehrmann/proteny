#!/usr/bin/python

import sys;
sys.path.append('utils');
sys.path.append('visualization');

import circos_data as cdata;
import circos_region as cr;
import circos_chr    as cc;
import proteny as ps;
import data;
import threshold as tr;

###############################################################################

savename = 'PR_schco_agabi_small.proteny';
circdir  = 'circos';

###############################################################################
# DATA

  # Download data for Schizophyllum commune
sc = data.example(data.schco2, 'scaffold_1');
  # Download data for Agaricus bisporus
ag = data.example(data.agabi, 'scaffold_2');

  # Prepare the protenty analysis
PR = ps.proteny();

  # Add our organisms.
id_a = PR.add_org(*sc, isfile=False);
id_b = PR.add_org(*ag, isfile=False);

###############################################################################
# ANALYSIS

k = PR.analyze(id_a=id_a, id_b=id_b);

  # Run BLAST
#k = PR.blast(id_a=id_a, id_b=id_b);
  # Calculate window scores for these results
#k = PR.hit_index(k);
  # Calculate genomic distances between ALL hits on different chromosomes
#k = PR.hit_distance(k);
  # Using a specified linkage, perform agglomerative clustering
#k = PR.hit_dendrogram(k);

  # Cut the dendrogram.
#k = PR.hit_cluster_height(k, 2000);
#k = PR.hit_cluster();

###############################################################################
# SAVE

  # Save our calculations
PR.save(savename);

###############################################################################
# VISUALIZE

      # name     [ (org_name, ID,       chrid,        start,   end,     strand) ... ]
  # Export data into a CIRCOS format
files = cdata.write_data(PR, k, '%s/data' % circdir);

  # Visualize a region
for reg in regions:
  cr.circos_region(files, reg, 30000, circdir, reg[0]);
#efor

  # Visualize relationships between chromosomes
agabichrs = ["agabi_%s" % s for s in PR.org_genomes[id_b].chrid() ];
for s in PR.org_genomes[id_a].chrid():
 cc.circos_chr(files, agabichrs + ["schco2_%s" % s], ["schco2_%s=0.4r" % s], circdir, "%s" % s );
#efor

