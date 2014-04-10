from ibidas import *

###############################################################################

# Each organism requires two files:
# Genes file:
#    |<chrid>\t<start>\t<end>\t<strand>\t<geneid>\t<transcriptid>\t<exonnumber>
# Sequence file:
#    |><chrid1>
#    |<sequence>
#    |><chrid2>
#    |<sequence>

__gene_slice_names__   = ( 'chrid', 'start', 'end', 'strand', 'geneid', 'transcriptid', 'exonid' );
__genome_slice_names__ = ( 'chrid', 'sequence');


###############################################################################
 
def agabi2(odir = None):
  genome_a = Read('Abisporus_varbisporusH97.v2.maskedAssembly.gz', format='fasta');
  genome_b = Read('Abisporus_var_bisporus.mitochondrion.scaffolds.fasta.gz', format='fasta');
  genome   = genome_a | Stack | genome_b;
  genes    = Read('Abisporus_varbisporusH97.v2.FilteredModels3.gff.gz');

  genes = genes[_.f2 == 'CDS'];
  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().split(' ')[1] for y in x.split(';')[1:]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1])).Detect();

  genes  = genes.To(_.f0, Do=_.Each(lambda x: x.split('_')[1]).Cast(str));
  genome = genome.To(_.f0, Do=_.Each(lambda x: x.split('_')[1]).Cast(str));

  genes  = genes / __gene_slice_names__;
  genome = genome / __genome_slice_names__;

  genes  = genes.Copy();
  genome = genome.Copy();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.fasta' % (odir));
  #fi

  return ("agabi2", genes, genome);
#edef

###############################################################################

def schco2(odir = None):
  genome = Read('Schco2_AssemblyScaffolds.fasta.gz', format='fasta');
  genes  = Read('Schco2_GeneCatalog_genes_20110923.gff.gz');

  genes = genes[_.f2 == 'CDS'];
  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().split(' ')[1] for y in x.split(';')[1:]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1])).Detect();

  genes  = genes.To(_.f0, Do=_.Each(lambda x: x.split('_')[1]).Cast(str));
  genome = genome.To(_.f0, Do=_.Each(lambda x: x.split('_')[1]).Cast(str));

  genes  = genes / __gene_slice_names__;
  genome = genome / __genome_slice_names__;

  genes  = genes.Copy();
  genome = genome.Copy();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.fasta' % (odir));
  #fi

  return ("schco2", genes, genome);
#edef

