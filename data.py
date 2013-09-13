dd = '/home/nfs/thiesgehrmann/groups/w/phd/data';

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
     

###############################################################################

org_names  = [];
org_genes  = [];
org_genome = [];

###############################################################################

def example(org, chr):
  (name, genes, genome) = org(odir = None);
  genes  = genes[_.f0 == chr];
  genome = genome[_.f0 == chr];

  return (name, genes, genome);
#edef

###############################################################################

def schco2(odir = '%s/1b-schco2' % dd):
  genome_a = Read(Fetch('http://genome.jgi-psf.org/Agabi_varbisH97_2/download/Abisporus_varbisporusH97.v2.maskedAssembly.gz'), format='fasta').Copy();
  genome_b = Read(Fetch('http://genome.jgi-psf.org/Agabi_varbisH97_2/download/Abisporus_var_bisporus.mitochondrion.scaffolds.fasta.gz')).Copy();
  genome   = genome_a | Stack | genome_b;
  genes    = Read(Fetch('http://genome.jgi-psf.org/Agabi_varbisH97_2/download/Abisporus_varbisporusH97.v2.FilteredModels3.gff.gz')).Copy()

  genes = genes[_.f2 == 'CDS'];
  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().split(' ')[1] for y in x.split(';')[1:]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1])).Detect();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.fasta' % (odir));
  #fi

  return ("schco2", genes, genome);
#edef

###############################################################################

def agabi(odir = '%s/2-agabi' % dd):
  genome = Read(Fetch('http://genome.jgi-psf.org/Schco2/download/Schco2_AssemblyScaffolds.fasta.gz')).Copy();
  genes  = Read(Fetch('http://genome.jgi-psf.org/Schco2/download/Schco2_GeneCatalog_genes_20110923.gff.gz')).Copy()
  
  genes = genes[_.f2 == 'CDS'];
  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().split(' ')[1] for y in x.split(';')[1:]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1])).Detect();


  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.fasta' % (odir));
  #fi

  return ("agabi", genes, genome);
#edef

###############################################################################

def human():
  odir = "%s/human" % dd;

  genome = Read('%s/Homo_sapiens.GRCh37.73.dna.fa' % (odir), sep=[' ']).Copy();
  genes  = Read('%s/Homo_sapiens.GRCh37.73.gtf' % (odir)).Copy();

  genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().replace('"', '').split(' ')[1] for y in x.split(';')[:-1]]));
  genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1]),_.f8.Each(lambda x: x[2]));

  genome.Get(_.f0, _.seq);

  Export(genes,  '%s/proteny_genes.tsv' % (odir));
  Export(genome, '%s/proteny_sequences.tsv' % (odir));

  return ("human", genes, genome);
#edef

###############################################################################

def mouse():
  odir = "%s/mus_musculus" % dd;

  genome = Read('%s/Mus_musculus.GRCm38.73.dna.fa' % (odir), sep=[' ']).Copy();
  genes  = Read('%s/Mus_musculus.GRCm38.73.gtf' % (odir)).Copy();

  genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().replace('"', '').split(' ')[1] for y in x.split(';')[:-1]]));
  genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1]),_.f8.Each(lambda x: x[2]));

  genome.Get(_.f0, _.seq);

  Export(genes,  '%s/proteny_genes.tsv' % (odir));
  Export(genome, '%s/proteny_sequences.tsv' % (odir));

  return ("mouse", genes, genome);
#edef

###############################################################################
