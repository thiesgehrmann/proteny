Proteny
=======

A tool to analyze synteny at the protein level.
We develop an algorithm to detect statistically significant clusters of exons between two proteomes.
The tool provides algorithms to quickly detect and visualize results in order to support conclusions from genomic data.

The method discovers clusters of hits from a bi-directional BLASTp of translated exon sequences in two organisms.
A dendrogram for hits is built based on genomic distances between hits, and cut based on significance of a cluster score based on a permutation test at each node in the tree.
The result is a set of large clusters describing high exonic conservation.

![An example of the figures generated by proteny](/readme/example_output.gif)

Algorithm
=========

We use BLASTp to produce a set of hits, which are used to build

![We use BLASTp to produce a set of hits, which are used to build](/readme/clustering_dendrogram_a.gif)

a dendrogram which is traversed to find

![a dendrogram which is traversed to find](/readme/clustering_dendrogram_b.gif)

significant clusters

![significant clusters.](/readme/clustering_dendrogram_c.gif)

Which are then visualized in a useful way

![Visualizations with Circos](/readme/visualization.gif)


Installation
=============

Proteny depends upon a few other packages.
(Circos is not necessary to run the program, but it is used to produce the visualizations)

 * Python: https://www.python.org/
 * IPython: http://ipython.org/
 * CIRCOS: http://circos.ca/
 * IBIDAS: https://github.com/mhulsman/ibidas
 * FastCluster: https://pypi.python.org/pypi/fastcluster
 * BLAST+: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/


1. Make sure that each dependency (and their dependencies) is installed.
2. Clone this github directory.
3. Modify the "$CIRCOS" variable in "circos_run" to reflect the location of your circos installation.

Running the example
=====================

An example is provided using data from the Yeast Gene Order Browser.
I produces Circos plots and outputs the clusters in a directory called 'example_output'.
It can be run with this line:

```shell
$> ./run_proteny.sh Cglabrata example_data/Cglabrata_genome.tsv example_data/Cglabrata_sequence.fasta Zrouxii example_data/Zrouxii_genome.tsv example_data/Zrouxii_sequence.fasta 0.05 1 example_output 10 2>&1 | tee proteny_out.log
```

Simple Usage
=============

The usage of the wrapper utility is given:

```shell
Run the proteny software
  Usage: ./run_proteny.sh <name_org1> <genes_org1> <genome_org1> <name_org2> <genes_org2> <genome_org2> <p-value> <c-thresh> <outdir>

  Arguments:
    <name_orgX>:   The short name for organism X
    <genes_orgX>:  The gene description file for organism X
    <genome_orgX>: The genome fasta file for organism X
    <p-value>:     The p-value threshold for clusters
    <c-thresh>:    The conservation threshold for clusters
    <outdir>:      The output directory to use
    <ncores>:      The number of cores available to the ipcluster
```

Input file specification
-------------------------

For each genome, you need to create two files:
 * A genes file, describing the exon locations of each gene
 * A fasta file, with the genome sequence of each chromosome.

The genes file must have the following columns, separated by columns:

 * **chromosome_id**: *string*, chromosome identifier
 * **start**: *integer*, the start location of this current exon
 * **end**: *integer*, the end location of this current exon
 * **strand**: *string*, the positve '+', or negative '-' strand
 * **gene_id**: *string*, a gene identifier
 * **transcript_id**: *string*, a transcript identifier
 * **exon_id**: *integer*, the exon number within the current transcript

The fasta file must have its sequence identifier corresponding to the chromosome id in the genes file.
Examples can be found in the example_data directory.


Advanced Usage
===============

Suppose you have these files for two organisms:
 * **ORG_NAME1**: genes_file_1.tsv, genome_file_1.fasta, and
 * **ORG_NAME2**: genes_file_2.tsv, genome_file_2.fasta,

then you can make a python file which runs Proteny in a more custimizable way:

```python

  # The null distribution generator for the significance test
from utils import cluster_null as null;

  # The visualization outputs
from visualization import circos_data as cdata;
from visualization import circos_chr  as cc;

  # The proteny functionality itself
import proteny as ps;
  # Functions to read the data
import data;

#######################################

  # A function that visualizes what we want,
  # in this case, all chromosomes for for org1 
  # are shown against those in org2, one by one.
def viz(PR, k, outdir):

  files = cdata.write_data(PR, k, '%s/data' % circdir);
  KS = PR.key_s(k);

  bchrs = ["%s_%s" % (PR.org_names[k['id_b']], str(i)) for i in PR.org_genomes[k['id_b']].Get(0)() ];
  for i in xrange(PR.org_genomes[k['id_a']].Shape()()):
    cc.circos_chr(files, bchrs + ["%s_%d" % (PR.org_names[k['id_a']], i+1)], ["%s_%d=0.4r" % (PR.org_names[k['id_a']], i+1)], circdir, "scaffold_%02d_%s" % (i+1, KS) );
    print "scaffold_%02d_%s" % (i+1, KS);
  #efor
#edef

#######################################

  # Read the data
org1 = data.read_prepared('ORG_NAME1', 'genes_file_1.tsv', 'genome_file_1.fasta');
org2 = data.read_prepared('ORG_NAME2', 'genes_file_2.tsv', 'genome_file_2.fasta');

  # Prepare the proteny structure
PR = ps.proteny();

  # Add our organisms.
id_a = PR.add_org(*org1, isfile=False);
id_b = PR.add_org(*org2, isfile=False);


  # Run analysis
k = PR.analyze(id_a=id_a, id_b=id_b,                     # Perform analysis between the two organisms we added
               cut='deeper_greater',                     # Find the smallest p-value given a conservation ratio
               nd=null.cluster_null_score_strict_smart,  # Our null distribution
               alpha=0.05,                               # p-value threshold
               ngenes_threshold=2,                       # Dont consider a cluster if it doesn't contain enough genes (not synteny)
               conservation_ratio=1);                    # Conservation ratio requirement

  # Save! It's important!
PR.save(savename);

  # Produce circos visualization files
viz(PR, k, circdir);

```
