dd = '/home/nfs/thiesgehrmann/groups/w/phd/data';

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

org_names =     [];
org_genes =     [];
org_sequences = [];

###############################################################################

def preprocess_data_schco2():
  odir = "%s/1b-schco2" % dd;

  sequences = Read('%s/schco2.assembly.fasta' % (odir));
  genes     = Read('%s/schco2.genes.gff' % (odir)).Copy()

  genes = genes[_.f2 == 'CDS'];
  genes = genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().split(' ')[1] for y in x.split(';')[1:]]));
  genes = genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1])).Detect();

  Export(genes,     '%s/proteny_genes.tsv' % (odir));
  Export(sequences, '%s/proteny_sequences.fasta' % (odir));
#edef

###############################################################################

def preprocess_data_agabi():
  odir = "%s/2-agabi" % dd;
  
  sequences = Read('%s/agabi.assembly_unmasked.fasta' % (odir));
  genes     = Read('%s/agabi.genes.gff' % (odir)).Copy()
  
  genes = genes[_.feature == 'CDS'];
  genes = genes.To(_.attribute, Do=_.Each(lambda x:[ y.strip().split(' ')[1] for y in x.split(';')[1:]]));
  genes = genes.Get(_.scaseqname, _.start, _.end, _.strand, _.attribute.Each(lambda x: x[0]), _.attribute.Each(lambda x: x[0]), _.attribute.Each(lambda x: x[1])).Detect();
  
  Export(genes,     '%s/proteny_genes.tsv' % (odir));
  Export(sequences, '%s/proteny_sequences.fasta' % (odir));
#edef


###############################################################################

def preprocess_data_human():
  odir = "%s/human" % dd;

  sequences = Read('%s/Homo_sapiens.GRCh37.73.dna.fa' % (odir), sep=[' ']).Copy();
  genes     = Read('%s/Homo_sapiens.GRCh37.73.gtf' % (odir)).Copy();

  genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().replace('"', '').split(' ')[1] for y in x.split(';')[:-1]]));
  genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1]),_.f8.Each(lambda x: x[2]));

  sequences.Get(_.f0, _.seq);

  Export(genes,     '%s/proteny_genes.tsv' % (odir));
  Export(sequences, '%s/proteny_sequences.tsv' % (odir));
#edef

###############################################################################

def preprocess_data_mouse():
  odir = "%s/mus_musculus" % dd;

  sequences = Read('%s/Mus_musculus.GRCm38.73.dna.fa' % (odir), sep=[' ']).Copy();
  genes     = Read('%s/Mus_musculus.GRCm38.73.gtf' % (odir)).Copy();

  genes.To(_.f8, Do=_.Each(lambda x:[ y.strip().replace('"', '').split(' ')[1] for y in x.split(';')[:-1]]));
  genes.Get(_.f0, _.f3, _.f4, _.f6, _.f8.Each(lambda x: x[0]), _.f8.Each(lambda x: x[1]),_.f8.Each(lambda x: x[2]));

  sequences.Get(_.f0, _.seq);

  Export(genes,     '%s/proteny_genes.tsv' % (odir));
  Export(sequences, '%s/proteny_sequences.tsv' % (odir));
#edef

###############################################################################
