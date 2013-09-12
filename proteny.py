import sys;
from ibidas import *;

sys.path.append('./utils');
import seq as sequtils;

###############################################################################

class proteny:

  ###############################################################################

  org_names   = [];
  org_genomes = [];
  org_genes   = [];
  org_exons   = [];

  blast_hits = {};

  ###############################################################################

  def __init__(self):
    self.org_names   = [];
    self.org_genomes = [];
    self.org_genes   = [];
    self.org_exons   = [];
    self.blast_hits  = {};
  #edef

  ###############################################################################

  def add_org(self, name, genes_filename, genome_filename, isrep=False):
    genome = self.load_genome(genome_filename, isrep=isrep);
    genes  = self.load_gene_defs(genes_filename, isrep=isrep);

    self.org_names.append(name);
    self.org_genomes.append(genome);
    self.org_genes.append(genes);
    self.org_exons.append(self.translate_exons(genome, genes));
  #edef
    

  ###############################################################################

  __gene_slice_names__ = ( 'chrid', 'start', 'end', 'strand', 'geneid', 'transcriptid', 'exonid');
  def load_gene_defs(self, data, isrep=False):
    return ((data if isrep else Read(data)) / self.__gene_slice_names__).Detect().Copy();
  #edef

  ###############################################################################

  __genome_slice_names__ = ('chrid', 'sequence');
  def load_genome(self, data, isrep=False):
    return ((data if isrep else Read(data)) / self.__genome_slice_names__).Detect().Copy();
  #edef

  #############################################################################

  __exon_slice_names__ = __gene_slice_names__ + ('sequence',)
  def translate_exons(self, genome, genes):

    EG = genes.GroupBy(_.transcriptid).Get(_.chrid[0], _.start, _.end, _.strand[0], _.geneid[0], _.transcriptid, _.exonid).Copy();

    dna = dict(zip(*genome()));

    exons = [];

    for gene in zip(*EG()):
      aid, starts, ends, strand, geneid, transcriptid, exonids = gene;
      chrdna = dna[aid];
      transcript = '';
      exonsource = [0];
      j = 0;
      for (s, e, id) in sorted(zip(starts, ends, exonids), key=lambda x: x[2], reverse=False):
        tr_ex       = chrdna[s-1:e];
        transcript += sequtils.revcomp(tr_ex) if (strand == '-') else tr_ex;
        exonsource += [exonsource[-1] + int((e-(s-1))/3)];
      #efor
      
      aa_seq = sequtils.translate(transcript);
      for (s, e, id) in sorted(zip(starts, ends, exonids), key=lambda x: x[2], reverse=False):
        ex_s = exonsource[id-1];
        ex_e = exonsource[id];
        exons.append((aid, s, e, strand, geneid, transcriptid, id, aa_seq[ex_s:ex_e]));
      #efor
    #efor

    return (Rep(exons) / self.__exon_slice_names__).Detect().Copy();
  #edef

  ###############################################################################

  def count_codons(self, genome, genes):

    EG = genes.GroupBy(_.transcriptid).Get(_.chrid[0], _.start, _.end, _.strand[0], _.geneid[0], _.transcriptid, _.exonid).Copy();

    dna = dict(zip(*genome()));

    cc = dict([(c, 0) for c in seq.codons]);
    total = 0;

    exons = [];

    for gene in zip(*EG()):
      aid, starts, ends, strand, geneid, transcriptid, exonids = gene;
      chrdna = dna[aid];
      transcript = '';
      exonsource = [0];
      j = 0;
      for (s, e, id) in sorted(zip(starts, ends, exonids), key=lambda x: x[2], reverse=False):
        tr_ex       = chrdna[s-1:e];
        transcript += sequtils.revcomp(tr_ex) if (strand == '-') else tr_ex;
        exonsource += [exonsource[-1] + int((e-(s-1))/3)];
      #efor
      for i in xrange(0, len(transcript), 3):
        codon = transcript[i:i+3].lower();
        if codon in cc:
          cc[codon] += 1;
          total     += 1;
        #fi
      #efor
    #efor

    return dict([(k, float(count) / float(total)) for (k,count) in cc.items() if len(k) == 3]);
  #edef

  ###############################################################################

  __blast_slice_names__ = ('a_chrid', 'a_strand', 'a_geneid', 'a_transcriptid', 'a_exonid', 'b_chrid', 'b_strand', 'b_geneid', 'b_transcriptid', 'b_exonid', 'qstart', 'qend', 'sstart', 'send', 'pident', 'mismatch', 'gapopen', 'evalue', 'bitscore');
  def blast(self, id_a = 0, id_b=1):
    if len(self.org_names) < 2:
      print "Cannot run BLAST, too few organisms added!";
    #fi

    a_exons = self.org_exons[id_a][_.end - _.start > 60] / tuple([ 'a_' + s for s in self.__exon_slice_names__]);
    b_exons = self.org_exons[id_b][_.end - _.start > 60] / tuple([ 'b_' + s for s in self.__exon_slice_names__]);
    
    print "Running BLAST for %s v %s" % (self.org_names[id_a], self.org_names[id_b]);
    R = a_exons | Blast(reciprocal = True, normalize=True, folder='./blast_runs/') | b_exons
    R = R.Copy();
    #filename = 'results.%s.%s.blast' % (self.org_names[id_a], self.org_names[id_b]);
    #self.blast_results.append((id_a, id_b, filename));
    #Save(R, filename);
    #R = Load(filename);
    
    F = R.Get(_.a_chrid, _.a_strand, _.a_geneid, _.a_transcriptid, _.a_exonid, \
              _.b_chrid, _.b_strand, _.b_geneid, _.b_transcriptid, _.b_exonid, \
              (_.a_start + (_.qstart * 3 )),                         \
              (_.a_start + (_.qend * 3 )),                           \
              (_.b_start + (_.sstart * 3 )),                         \
              (_.b_start + (_.send * 3 )),                           \
              _.pident,                                              \
              _.mismatch,                                            \
              _.gapopen,                                             \
              _.evalue,                                              \
              _.bitscore);
    F = ( F / self.__blast_slice_names__).Copy()
    self.blast_hits[(id_a, id_b)] = F;
  #edef

  ###############################################################################

  
  

###############################################################################


