import sys;
import types;
import cPickle;

from ibidas import *;
from scipy.cluster import hierarchy;
import numpy as np;

sys.path.append('./utils');
import seq as sequtils;
import smoothing;

###############################################################################

class proteny:

  #############################################################################

  org_names   = [];
  org_genomes = [];
  org_genes   = [];
  org_exons   = [];

  blast_hits       = {};
  hit_windows      = {};
  hit_distances    = {};
  cluster_linkages = {};
  hit_clusters     = {};

  ###############################################################################

  def __init__(self, load=None):
    self.org_names   = [];
    self.org_genomes = [];
    self.org_genes   = [];
    self.org_exons   = [];

    self.blast_hits       = {};
    self.hit_windows      = {};
    self.hit_distances    = {};
    self.cluster_linkages = {};
    self.hit_clusters     = {};

    if not(load == None):
      p = cPickle.load(open(load, 'r'));
      for (attr, value) in p:
        setattr(self, attr, value);
      #efor
    #fi
  #edef

  ###############################################################################

  def save(self, file):
    fd = open(file, 'w');
    isfunc = lambda obj, attr: hasattr(obj, attr) and type(getattr(obj, attr)) == types.MethodType;
    p = [ (attr, getattr(self, attr)) for attr in dir(self) if not(isfunc(self, attr)) ];

    cPickle.dump(p, fd, protocol=2);
    fd.close();
  #edef

  #############################################################################

  def add_org(self, name, genes_data, genome_data, isfile=True):
    genome = self.load_genome(genome_data, isfile=isfile);
    genes  = self.load_gene_defs(genes_data, isfile=isfile);

    self.org_names.append(name);
    self.org_genomes.append(genome);
    self.org_genes.append(genes);
    self.org_exons.append(self.translate_exons(genome, genes));

    return len(self.org_names)-1;
  #edef

  ###############################################################################

  __gene_slice_names__ = ( 'chrid', 'start', 'end', 'strand', 'geneid', 'transcriptid', 'exonid');

  def load_gene_defs(self, data, isfile=True):
    return ((Read(data) if isfile else data) / self.__gene_slice_names__).Detect().Copy();
  #edef

  ###############################################################################

  __genome_slice_names__ = ('chrid', 'sequence');

  def load_genome(self, data, isfile=True):
    return ((Read(data) if isfile else data) / self.__genome_slice_names__).Detect().Copy();
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

  __blast_slice_names__ = ('a_chrid', 'a_strand', 'a_geneid', 'a_transcriptid', 'a_exonid', \
                           'b_chrid', 'b_strand', 'b_geneid', 'b_transcriptid', 'b_exonid', \
                           'a_start', 'a_end', 'b_start', 'b_end',                          \
                           'pident', 'mismatch', 'gapopen', 'evalue', 'bitscore');

  def blast(self, id_a = 0, id_b=1):
    if len(self.org_names) < 2:
      print "Cannot run BLAST, too few organisms added!";
    #fi

    a_exons = self.org_exons[id_a][_.end - _.start > 60] / tuple([ 'a_' + s for s in self.__exon_slice_names__]);
    b_exons = self.org_exons[id_b][_.end - _.start > 60] / tuple([ 'b_' + s for s in self.__exon_slice_names__]);
    
    print "Running BLAST for %s v %s" % (self.org_names[id_a], self.org_names[id_b]);
    R = a_exons | Blast(normalize=True, folder='./blast_runs/') | b_exons
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

    return (id_a, id_b);
  #edef

  ###############################################################################

  def windows(self, k, WS=[20000]):

    BR = self.blast_hits[k];
    F  = smoothing.reindex_blast(BR);

    hits = F.Get(_.i, _.a_chrid, _.a_geneid, _.a_exonid, _.b_chrid, _.b_geneid, _.b_exonid, _.pident, _.evalue, _.bitscore);
    BR_a = F.Get(_.i, _.a_chrid, _.a_start, _.a_end) / ('i', 'chrid', 'start', 'end');
    BR_b = F.Get(_.i, _.b_chrid, _.b_start, _.b_end) / ('i', 'chrid', 'start', 'end');

    brs  = [ BR_a, BR_b ];

    H = [ [[]] + list(x[1:]) for x in zip(*hits()) ]
    O = [];

    for br in brs:
      org_ind, H = smoothing.sort_org(br, H);
      O.append(org_ind);
    #efor

    RS = [];
    hlen = len(H);

    print "Smoothing hits. This will take a while!";
    for ws in WS:
      R = [];
      for h in xrange(hlen):
        print "%d/%d" % (h, hlen);
        R.append(smoothing.get_reg_hits(H, O, h, ws));
      #efor
      RS.append(R);
    #efor

    scores = [[ smoothing.score(RS[i][j]) for i in xrange(len(WS)) ] for j in xrange(hlen)]

    self.hit_windows[k] = (H, O, scores, WS);

    return k;
  #edef

  #############################################################################

  def hit_distance(self, k):
    H, O, scores, WS = self.hit_windows[k];

    HC = smoothing.chr_pair_group(H);
    D  = {};

    i = 0;
    for hc_k in HC.keys():
      print i, len(HC.keys());
      hc  = HC[hc_k];
      if len(hc) < 2:
        continue;
      #fi
      cdm = np.array(smoothing.condenseddm(H, O, hc), dtype=np.dtype('u8'));
      if len(cdm) > 2:
        D[hc_k] = cdm;
      #fi
      i += 1;
    #efor
    self.hit_distances[k] = (HC, D);
    return k;
  #edef

  #############################################################################

  def cluster_linkage(self, k, linkage_type='single'):
    if k not in self.hit_distances:
      print "You must run hit_distance() first!";
      return None;
    #fi

    linkage_types = { 'single'   : hierarchy.single,
                      'complete' : hierarchy.complete,
                      'average'  : hierarchy.average,
                      'weighted' : hierarchy.weighted,
                      'centroid' : hierarchy.centroid,
                      'median'   : hierarchy.median,
                      'ward'     : hierarchy.ward };

    HC, D = self.hit_distances[k];
    L = {};
    for dk in D.keys():
      print dk;
      L[dk] = linkage_types[linkage_type](D[dk]);
    #efor

    self.cluster_linkages[k] = (L, linkage_type);

    return k;
  #edef

  #############################################################################

  def cluster_hits(self, k, CS=[20000]):
    if (k not in self.hit_windows) or (k not in self.cluster_linkages):
      print "You must run windows() and cluster_linkage() first!";
      return None;
    #fi

    H,  O, scores, WS = self.hit_windows[k];
    HC, D             = self.hit_distances[k];
    L,  lt            = self.cluster_linkages[k];

    CD = [];
    for cs in CS:

      HC_clust = {};
      for lk in L.keys():
        linkage = L[lk];
        hc      = HC[lk];

        hclust  = hierarchy.fcluster(linkage, cs, criterion='distance');
        nclust  = max(hclust);
        clusts  = [ [] for i in xrange(nclust) ];
        
        for (h, c) in zip(hc, hclust):
          clusts[c-1].append(h);
        #efor
        HC_clust[lk] = clusts;
      #efor

      hit_clusts = [];
      for lk in HC_clust.keys():
        clusts = HC_clust[lk];
        for C in clusts:
          cd = smoothing.clust_description(H, O, scores, C);
          hit_clusts.append(cd);
        #efor
      #efor
      CD.append(hit_clusts);
    #efor

    self.hit_clusters[k] = (CD, CS);
    return k;
  #edef

  #############################################################################
  

#eclass

###############################################################################

