#!/bin/env python

import sys;
import os;

import cluster_null as null;
#from utils import cluster_null as null;

from visualization import circos_data as cdata;
from visualization import circos_region as cr;
from visualization import circos_chr    as cc;

import data;
import proteny as ps;
import postprocessing as pp;

sys.path= sys.path + [ os.path.dirname(os.path.realpath(__file__)) + ('/%s/' % x if x != '' else '')  for x in ['visualization'] ];

###############################################################################

def usage(arg0):
  print "Run the core Protreny algorithm";
  print "Usage %s <name_org1> <genes_org1> <genome_org1> <name_org2> <genes_org2> <genome_org2> <p-value> <c-thresh> <outdir>";
  print "";
  print "  Arguments:";
  print "    <name_orgX>:   The short name for organism X";
  print "    <genes_orgX>:  The gene description file for organism X";
  print "    <genome_orgX>: The genome fasta file for organism X";
  print "    <p-value>:     The p-value threshold for clusters";
  print "    <c-thresh>:    The conservation threshold for clusters";
  print "    <outdir>:      The output directory to use";
 
#edef

###############################################################################

def viz(PR, k, outdir):
  KS = PR.key_s(k);

  files = cdata.write_data(PR, k, '%s/data' % outdir);

  # Visualize relationships between chromosomes
  bchrs = ["%s_%s" % (PR.org_names[k['id_b']], str(i)) for i in PR.org_genomes[k['id_b']].Get(0)() ];
  for i in xrange(PR.org_genomes[k['id_a']].Shape()()):
    cc.circos_chr(files, bchrs + ["%s_%d" % (PR.org_names[k['id_a']], i+1)], ["%s_%d=0.4r" % (PR.org_names[k['id_a']], i+1)], outdir, "scaffold_%02d_%s" % (i+1, KS) );
    print "scaffold_%02d_%s" % (i+1, KS);
  #efor

#edef

###############################################################################

if __name__ == '__main__':

  if len(sys.argv) < 10:
    usage(sys,argv[0]);
    sys.exit(1);
  #fi

  name_org1   = sys.argv[1];
  genes_org1  = sys.argv[2];
  genome_org1 = sys.argv[3];

  name_org2   = sys.argv[4];
  genes_org2  = sys.argv[5];
  genome_org2 = sys.argv[6];
  
  pvalue  = float(sys.argv[7]);
  cthresh = float(sys.argv[8]);
  outdir  = sys.argv[9];

  savename = '%s/proteny_data.proteny' % outdir;

    # Load the data
  org1 = data.read_prepared(name_org1, genes_org1, genome_org1);
  org2 = data.read_prepared(name_org2, genes_org2, genome_org2);

    # Prepare proteny structure
  PR = ps.proteny();

    # Add our organisms.
  id_a = PR.add_org(*org1, isfile=False);
  id_b = PR.add_org(*org2, isfile=False);
  
  ###############################################################################
  # ANALYSIS

    # Run analysis
  k = PR.analyze(id_a=id_a, id_b=id_b, cut='deeper_greater', nd=null.cluster_null_score_strict_smart, alpha=pvalue, ngenes_threshold=2, conservation_ratio=cthresh);

    # Save! It's important!
  PR.save(savename);

    # Produce circos visualizations
  viz(PR, k, outdir);

#fi
