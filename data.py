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

def example(org, chr, odir=None):
  (name, genes, genome) = org(odir = odir);
  genes  = genes[_.f0 == chr];
  genome = genome[_.f0 == chr];

  return (name, genes.Copy(), genome.Copy());
#edef

###############################################################################

def agabi(odir = None):
  genome_a = Read(Fetch('http://genome.jgi-psf.org/Agabi_varbisH97_2/download/Abisporus_varbisporusH97.v2.maskedAssembly.gz'), format='fasta');
  genome_b = Read(Fetch('http://genome.jgi-psf.org/Agabi_varbisH97_2/download/Abisporus_var_bisporus.mitochondrion.scaffolds.fasta.gz'));
  genome   = genome_a | Stack | genome_b;
  genes    = Read(Fetch('http://genome.jgi-psf.org/Agabi_varbisH97_2/download/Abisporus_varbisporusH97.v2.FilteredModels3.gff.gz'));

  genes = genes[_.f2 == 'CDS'];
  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().split(' ')[1] for y in x.split(';')[1:]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1])).Detect();

  genes  = genes.Copy();
  genome = genome.Copy();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.fasta' % (odir));
  #fi

  return ("agabi", genes, genome);
#edef

###############################################################################

def schco2(odir = None):
  genome = Read(Fetch('http://genome.jgi-psf.org/Schco2/download/Schco2_AssemblyScaffolds.fasta.gz'));
  genes  = Read(Fetch('http://genome.jgi-psf.org/Schco2/download/Schco2_GeneCatalog_genes_20110923.gff.gz'));
  
  genes = genes[_.f2 == 'CDS'];
  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().split(' ')[1] for y in x.split(';')[1:]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1])).Detect();

  genes  = genes.Copy();
  genome = genome.Copy();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.fasta' % (odir));
  #fi

  return ("schco2", genes, genome);
#edef

###############################################################################

def human(odir = "%s/human" % dd):

  genome = Read('%s/human/Homo_sapiens.GRCh37.73.dna.fa' % (dd), sep=[' ']);
  genes  = Read('%s/human/Homo_sapiens.GRCh37.73.gtf'    % (dd));

  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().replace('"', '').split(' ')[1] for y in x.split(';')[:-1]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1]),_.f8.Each(lambda x: x[2]));

  genome = genome.Get(_.f0, _.seq);

  genes  = genes.Copy();
  genome = genome.Copy();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.tsv' % (odir));
  #fi

  return ("human", genes, genome);
#edef

###############################################################################

def mouse(odir = "%s/mus_musculus" % dd):
  genome = Read('%s/mus_musculus/Mus_musculus.GRCm38.73.dna.fa' % (dd), sep=[' ']);
  genes  = Read('%s/mus_musculus/Mus_musculus.GRCm38.73.gtf'    % (dd));

  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().replace('"', '').split(' ')[1] for y in x.split(';')[:-1]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1]),_.f8.Each(lambda x: x[2]));

  genome = genome.Get(_.f0, _.seq);

  genes  = genes.Copy();
  genome = genome.Copy();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.tsv' % (odir));
  #fi

  return ("mouse", genes, genome);
#edef

###############################################################################

def aniger_n402(odir = None):
  genome = Read('/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/n402_sequence/assembly/n402_atcc.unpadded.fasta', sep=[]);
  genes  = Read('/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/n402_sequence/annotations/results/n402_annotations.gff');

  #<chrid>\t<start>\t<end>\t<strand>\t<geneid>\t<transcriptid>\t<exonnumber>
  genes = genes[_.f2 == 'exon'];
  genes = genes.To(_.f8, Do=_.Each(lambda x: x.split(';')));
  genes = genes.To(_.f0, Do=_.Each(lambda x: x.split('_')[0]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[1].split('=')[1].split('_')[0]), _.f8.Each(lambda x: int(x[0].split('=')[1].split('_')[1][4:])+1 ));
  genes = genes / ("f0", "f1", "f2", "f3", "f4", "f5");
  genes = genes.Get(_.f0, _.f1, _.f2, _.f3, _.f4, _.f4, _.f5);

  genome = genome.To(_.f0, Do=_.Each(lambda x: x.split('_')[0]));

  genes  = genes.Detect().Copy();
  genome = genome.Detect().Copy();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.tsv' % (odir));
  #fi

  return ("aniger_n402", genes, genome);
#edef

###############################################################################

def aniger_513_88(odir = None):
  genome = Read(Fetch('http://www.aspergillusgenome.org/download/sequence/A_niger_CBS_513_88/current/A_niger_CBS_513_88_current_chromosomes.fasta.gz'));
  genes  = Read(Fetch('http://www.aspergillusgenome.org/download/gff/A_niger_CBS_513_88/A_niger_CBS_513_88_current_features.gff'));

  genes = genes[_.f2 == 'exon'];
  genes = genes.To(_.f8, Do=_.Each(lambda x: x.split(';')));
  genes = genes.To(_.f0, Do=_.Each(lambda x: x.split('_')[0]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[1].split('=')[1].split('-')[0]), _.f8.Each(lambda x: x[0].split('-')[2][1:]))
  genes = genes / ("f0", "f1", "f2", "f3", "f4", "f5");
  genes = genes.Get(_.f0, _.f1, _.f2, _.f3, _.f4, _.f4, _.f5);

  genome = genome.To(_.f0, Do=_.Each(lambda x: x.split('_')[0]));

  genes  = genes.Detect().Copy();
  genome = genome.Detect().Copy();

  if not(odir == None):
    Export(genes,  '%s/proteny_genes.tsv' % (odir));
    Export(genome, '%s/proteny_sequences.tsv' % (odir));
  #fi

  return ("aniger_513_88", genes, genome);
#edef

###############################################################################
