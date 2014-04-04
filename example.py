#!/usr/bin/python

import sys;
import os;

from utils import cluster_null as null;

from visualization import circos_data as cdata;
from visualization import circos_region as cr;
from visualization import circos_chr    as cc;

import proteny as ps;
import data;
import postprocessing as pp;

###############################################################################
# VISUALIZE

def viz(PR, k, circdir):
  KS = PR.key_s(k);

        # name     [ (org_name, ID,       chrid,        start,   end,     strand) ... ]
  regions = \
       [('bri1',   [ ('schco2', '1',  5346259, 5348750), ('agabi2', '1',  1865954, 1869064) ] ),
        ('c2h2',   [ ('schco2', '16', 363313,  364601 ), ('agabi2', '18', 156741,  158622 ) ] ),
        ('fst3',   [ ('schco2', '6',  1486211, 1492035), ('agabi2', '13', 776501,  779982 ) ] ),
        ('fst4',   [ ('schco2', '4',  2515305, 2519027), ('agabi2', '7',  344969,  349128 ) ] ),
        ('gat1',   [ ('schco2', '1',  1406704, 1408481), ('agabi2', '14', 113077,  114790 ) ] ),
        ('hom1',   [ ('schco2', '7',  1749469, 1751039), ('agabi2', '5',  549916,  551732 ) ] ),
        ('hom2',   [ ('schco2', '8',  1147947, 1150388), ('agabi2', '5',  1523421, 1524802) ] ),
        ('wc1',    [ ('schco2', '9',  1046213, 1049419), ('agabi2', '8',  25945,   28694  ) ] ),
        ('wc2',    [ ('schco2', '2',  844864,  846181 ), ('agabi2', '2',  875026,  876315 ) ] ),
        ('random', [ ('schco2', '5',  2483079, 2483664), ('agabi2', '8',  177699,  179015 ) ] ) ];
  interesting_genes = [ ( 'bri1', 255701 ), 
                        ( 'c2h2', 114363 ),
                        ( 'fst3', 257422),
                        ( 'fst4', 66861),
                        ( 'gat1', 255004),
                        ( 'hom1', 257652),
                        ( 'hom2', 1034289),
                        ( 'wc1',  78657),
                        ( 'wc2',  13988) ]

    # Export data into a CIRCOS format
  files = cdata.write_data(PR, k, '%s/data' % circdir);

    # Clusters containing our genes.
  fc = pp.hit_clusters_containing(PR, k, [ j[1] for j in interesting_genes ], [])[0];
  for (g, rs) in zip([ j[0] for j in interesting_genes], fc):
    for r in rs:
      reg = pp.hit_clust_2_reg(PR, k, r, name='%s_%d' % (g, r));
      cr.circos_region(files, reg, 100, circdir, 'reg_discovered_' + KS + '_' + reg[0]);
      print 'RegDiscovered_' + KS + '_' + reg[0] + '.conf';
    #efor
  #efor

    # Overlapping clusters
  IR = pp.hit_cluster_overlap(PR, k);
  IR = sorted(IR, key=lambda x: len(x[1]));
  IRs = [ ir for ir in IR if ir[0] == '2' ];

  for reg in IRs:
    cr.circos_region(files, reg, 100, circdir, ('RegOverlapping_' + KS + '_' + reg[0]));
    print 'RegOverlapping_' + KS + '_' + reg[0];
  #efor


    # Visualize a region
  for reg in regions:
    cr.circos_region(files, reg, 30000, circdir, ('RegKnown_' + KS + '_' + reg[0]));
    print 'RegKnown_' + KS + '_' + reg[0];
  #efor

  # Visualize relationships between chromosomes
  bchrs = ["%s_%s" % (PR.org_names[k['id_b']], str(i)) for i in PR.org_genomes[k['id_b']].Get(0)() ];
  for i in xrange(PR.org_genomes[k['id_a']].Shape()()):
    cc.circos_chr(files, bchrs + ["%s_%d" % (PR.org_names[k['id_a']], i+1)], ["%s_%d=0.4r" % (PR.org_names[k['id_a']], i+1)], circdir, "scaffold_%02d_%s" % (i+1, KS) );
    print "scaffold_%02d_%s" % (i+1, KS);
  #efor


#  bchrs = ["agabi2_%d" % (i+1) for i in xrange(PR.org_genomes[k['id_b']].Shape()()) ];
#  for i in xrange(36):
#    cc.circos_chr(files, bchrs + ["schco2_%d" % (i+1)], ["schco2_%d=0.4r" % (i+1)], circdir, "scaffold_%02d_%s" % (i+1, KS) );
#    print "scaffold_%02d_%s" % (i+1, KS);
#  #efor

#edef

###############################################################################

savename = 'PROTENY_OUTPUT_basid.proteny';
circdir  = 'PROTENY_OUTPUT/basid';

###############################################################################
# DATA

  # Download data for Schizophyllum commune
sc = data.schco2();
  # Download data for Agaricus bisporus
ag = data.agabi2();

  # Prepare the protenty analysis
PR = ps.proteny();

  # Add our organisms.
id_a = PR.add_org(*sc, isfile=False);
id_b = PR.add_org(*ag, isfile=False);

###############################################################################
# ANALYSIS

  # Run analysis
k = PR.analyze(id_a=id_a, id_b=id_b, cut='deeper_greater', nd=null.cluster_null_score_strict_smart, alpha=0.05, ngenes_threshold=2, conservation_ratio=1);
  # Save! It's important!
PR.save(savename);

  # Produce circos visualizations
viz(PR, k, circdir);

#for r in [ 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0 ]:
#  k = PR.hit_cluster(k, cut='deeper_greater', nd=null.cluster_null_score_strict_smart, conservation_ratio=r)
#  fout = 'PR_basid_r_test_%0.2f' % r;
#  circdir = 'visualizations/basid_moreperms_check/complete_deeper_greater_%0.2f' % r;
#  PR.save(fout);
#  viz(PR, k);
#efor

