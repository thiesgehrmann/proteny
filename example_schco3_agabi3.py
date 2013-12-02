#!/usr/bin/python

import sys;
sys.path.append('utils');
sys.path.append('visualization');

import cluster_null as null;
import circos_data as cdata;
import circos_region as cr;
import circos_chr    as cc;
import proteny as ps;
import data;

###############################################################################

savename = 'PR_schco3_agabi3.proteny';
circdir  = 'circos';

###############################################################################
# DATA

  # Download data for Schizophyllum commune
sc = data.schco3();
  # Download data for Agaricus bisporus
ag = data.agabi3();

  # Prepare the protenty analysis
PR = ps.proteny();

  # Add our organisms.
id_a = PR.add_org(*sc, isfile=False);
id_b = PR.add_org(*ag, isfile=False);

###############################################################################
# VISUALIZE

def viz(PR, k):
  KS = PR.key_s(k);

        # name     [ (org_name, ID,       chrid,        start,   end,     strand) ... ]
  regions = \
       [('bri1',   [ ('schco2', 'scaffold_1',  5346259, 5348750), ('agabi', 'scaffold_1',  1865954, 1869064) ] ),
        ('c2h2',   [ ('schco2', 'scaffold_16', 363313,  364601 ), ('agabi', 'scaffold_18', 156741,  158622 ) ] ),
        ('fst3',   [ ('schco2', 'scaffold_6',  1486211, 1492035), ('agabi', 'scaffold_13', 776501,  779982 ) ] ),
        ('fst4',   [ ('schco2', 'scaffold_4',  2515305, 2519027), ('agabi', 'scaffold_7',  344969,  349128 ) ] ),
        ('gat1',   [ ('schco2', 'scaffold_1',  1406704, 1408481), ('agabi', 'scaffold_14', 113077,  114790 ) ] ),
        ('hom1',   [ ('schco2', 'scaffold_7',  1749469, 1751039), ('agabi', 'scaffold_5',  549916,  551732 ) ] ),
        ('hom2',   [ ('schco2', 'scaffold_8',  1147947, 1150388), ('agabi', 'scaffold_5',  1523421, 1524802) ] ),
        ('wc1',    [ ('schco2', 'scaffold_9',  1046213, 1049419), ('agabi', 'scaffold_8',  25945,   28694  ) ] ),
        ('wc2',    [ ('schco2', 'scaffold_2',  844864,  846181 ), ('agabi', 'scaffold_2',  875026,  876315 ) ] ),
        ('random', [ ('schco2', 'scaffold_5',  2483079, 2483664), ('agabi', 'scaffold_8',  177699,  179015 ) ] ) ];
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
      cr.circos_region(files, reg, 0, circdir, 'reg_discovered_' + KS + '_' + reg[0]);
      print 'RegDiscovered_' + KS + '_' + reg[0] + '.conf';
    #efor
  #efor

    # Overlapping clusters
  IR = pp.hit_cluster_overlap(PR, k);
  IR = sorted(IR, key=lambda x: len(x[1]));
  IR = [ IR[-1] ];

  for reg in IR:
    cr.circos_region(files, reg, 100, circdir, ('RegOverlapping_' + KS + '_' + reg[0]));
    print 'RegOverlapping_' + KS + '_' + reg[0];
  #efor


    # Visualize a region
  for reg in regions:
    cr.circos_region(files, reg, 30000, circdir, ('RegKnown_' + KS + '_' + reg[0]));
    print 'RegKnown_' + KS + '_' + reg[0];
  #efor

  # Visualize relationships between chromosomes
  agabichrs = ["agabi_scaffold_%d" % (i+1) for i in xrange(21) ];
  for i in xrange(36):
    cc.circos_chr(files, agabichrs + ["schco2_scaffold_%d" % (i+1)], ["schco2_scaffold_%d=0.4r" % (i+1)], circdir, "scaffold_%02d_%s" % (i+1, KS) );
    print "Scaffold_%02d_%s" % (i+1, KS);
  #efor

#edef

###############################################################################
# ANALYSIS

k = PR.analyze(id_a=id_a, id_b=id_b);

nulls = [ null.cluster_null_score_gentle ]; 
cuts  = [ 'simple' ];

for N in nulls:
  for C in cuts:
    k = PR.hit_cluster(k, alpha=0.05, cut=C, nd=N);
    PR.save(savename);
    print k;
    #viz(PR, k);
  #efor
#efor    

###############################################################################
# SAVE

  # Save our calculations
PR.save(savename);

