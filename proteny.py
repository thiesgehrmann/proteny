import sys;
import types;
import cPickle;

from ibidas import *;
from scipy.cluster import hierarchy;
import numpy as np;

sys.path.append('./utils');
import seq as sequtils;
import smoothing;
import cluster as cluster;
reload(cluster);
import util as util;
reload(util);

from ibidas.utils.util import debug_here;

###############################################################################

class proteny:

  #############################################################################

  org_names   = [];
  org_genomes = [];
  org_genes   = [];
  org_exons   = [];
  org_chrs    = [];

  blast_hits          = {};
  hits                = {};
  hit_window_index    = None;
  hit_windows         = {};
  hit_distances       = {};
  hit_dendrograms     = {};
  hit_clusters        = {};
  hit_clusters_height = {};

  ###############################################################################

  def __init__(self, load=None):
    self.org_names   = [];
    self.org_genomes = [];
    self.org_genes   = [];
    self.org_exons   = [];
    self.org_chrs    = [];

    self.blast_hits          = {};
    self.hits                = {};
    self.hit_window_index    = None;
    self.hit_windows         = {};
    self.hit_distances       = {};
    self.hit_dendrograms     = {};
    self.hit_clusters        = {};
    self.hit_clusters_height = {};

    if not(load == None):
      p = cPickle.load(open(load, 'r'));
      for (attr, value) in p:
        setattr(self, attr, value);
      #efor
    #fi
  #edef

  ###############################################################################

  def analyze(self, id_a=0, id_b=1):

      # Run BLAST
    k = self.blast(id_a=id_a, id_b=id_b);
      # Calculate window scores for these results
    k = self.hit_index(k);
      # Calculate genomic distances between ALL hits on different chromosomes
    k = self.hit_distance(k);
      # Using a specified linkage, perform agglomerative clustering
    k = self.hit_dendrogram(k);

      # Cut the dendrogram.
    k = self.hit_cluster(k);

    return k;
  #edef

  ###############################################################################

  def key(self, k):
    key_elems = [ 'id_a' , 'id_b', 'linkage_type', 'H', 'alpha' ];
    nk        = dict([ (f, None) for f in key_elems]);

    for (i,v) in enumerate(k):
      nk[key_elems[i]] = v;
    #efor
    return nk;
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
    id = len(self.org_names);

    genome = self.load_genome(genome_data, isfile=isfile);
    genes  = self.load_gene_defs(genes_data, isfile=isfile);

    self.org_names.append(name);
    self.org_genomes.append(genome);
    self.org_genes.append(genes);
    self.org_exons.append(self.translate_exons(genome, genes));
    self.org_chrs.append(util.prep_exon_list(self.org_exons[id]));

    return id;
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

    cc = dict([(c, 0) for c in sequtils.codons]);
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
                           'coverage', 'score', 'pident', 'mismatch', 'gapopen', 'evalue', 'bitscore');

  def blast(self, id_a = 0, id_b=1):
    if len(self.org_names) < 2:
      print "Cannot run BLAST, too few organisms added!";
      return None;
    #fi

    a_exons = self.org_exons[id_a][_.end - _.start > 60] / tuple([ 'a_' + s for s in self.__exon_slice_names__]);
    a_exons = a_exons.To(_.a_sequence, Do=_.ReplaceMissing());
    b_exons = self.org_exons[id_b][_.end - _.start > 60] / tuple([ 'b_' + s for s in self.__exon_slice_names__]);
    b_exons = b_exons.To(_.b_sequence, Do=_.ReplaceMissing());
    
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
              ((_.qend.Cast(float) - _.qstart + 1) * 3 + (_.send.Cast(float) - _.sstart + 1) * 3) / ((_.qlen.Cast(float) + _.slen) * 3), \
              ((_.qend.Cast(float) - _.qstart + 1) * 3 + (_.send.Cast(float) - _.sstart + 1) * 3) / ((_.qlen.Cast(float) + _.slen) * 3) * _.evalue.Each(lambda x: 1.0 - min(1.0,x)), \
              _.pident,                                              \
              _.mismatch,                                            \
              _.gapopen,                                             \
              _.evalue,                                              \
              _.bitscore);
    F = ( F / self.__blast_slice_names__).Copy()
    self.blast_hits[(id_a, id_b)] = F;

    return self.key((id_a, id_b));
  #edef

  ###############################################################################

  def windows(self, k, WS=0):

    k['WS'] = WS;

    if self.hit_window_index != None:
      H, O = self.hit_window_index;
    else:
      if (k['id_a'], k['id_b']) not in self.blast_hits:
        print "You need to run blast() first!";
        return None;
      #fi

      BR = self.blast_hits[(k['id_a'], k['id_b'])];
      F  = util.reindex_blast(BR);

      hits = F.Get(_.i, _.a_chrid, _.a_geneid, _.a_exonid, _.b_chrid, _.b_geneid, _.b_exonid, _.coverage, _.pident, _.evalue, _.bitscore);
      BR_a = F.Get(_.i, _.a_chrid, _.a_start, _.a_end) / ('i', 'chrid', 'start', 'end');
      BR_b = F.Get(_.i, _.b_chrid, _.b_start, _.b_end) / ('i', 'chrid', 'start', 'end');

      brs  = [ BR_a, BR_b ];

      H = [ [[]] + list(x[1:]) for x in zip(*hits()) ]
      O = [];

      for br in brs:
        org_ind, H = smoothing.sort_org(br, H);
        O.append(org_ind);
      #efor
      self.hit_window_index = (H, O);
    #fi

    hlen = len(H);

    print "Smoothing hits. This will take a while!";
    RS = [];
    for h in xrange(hlen):
      print "\r%d/%d" % (h+1, hlen),
      sys.stdout.flush();
      RS.append(smoothing.get_reg_hits(H, O, h, WS));
    #efor
    print "";

    scores = [ smoothing.score(RS[j]) for j in xrange(hlen)]

    self.hit_windows[(k['id_a'], k['id_b'], k['WS'])] = scores;

    return k;
  #edef

  #############################################################################

  def hit_index(self, k):
    BR = self.blast_hits[(k['id_a'], k['id_b'])];
    F  = util.reindex_blast(BR);

    self.hits[(k['id_a'], k['id_b'])] = zip(*F());

    return k;
  #edef

  #############################################################################

  def hit_distance(self, k):
    if (k['id_a'], k['id_b']) not in self.hits:
      print "You must run hit_index() first!";
      return None;
    #fi

    hits = self.hits[(k['id_a'], k['id_b'])];

    HC, D = cluster.calc_distances(hits);

    self.hit_distances[(k['id_a'], k['id_b'])] = (HC, D);
    
    return k;
  #edef

  #############################################################################

  def hit_dendrogram(self, k, linkage_type='single'):
    k['linkage_type'] = linkage_type;
    
    if (k['id_a'], k['id_b']) not in self.hit_distances:
      print "You must run hit_distance() first!";
      return None;
    #fi
    
    HC, D = self.hit_distances[(k['id_a'], k['id_b'])];

    T = cluster.calc_dendrograms(HC, D, linkage_type);
    
    self.hit_dendrograms[(k['id_a'], k['id_b'], k['linkage_type'])] = T;
    
    return k;
  #edef

  #############################################################################

  def hit_cluster(self, k, alpha=0.05):

    if ((k['id_a'], k['id_b']) not in self.hits) or \
       ((k['id_a'], k['id_b']) not in self.hit_distances) or \
       ((k['id_a'], k['id_b'], k['linkage_type']) not in self.hit_dendrograms):
      print "You must run hit_index(), hit_distance() and hit_dendrogram() first!";
      return None;
    #fi

    k['alpha'] = alpha;

    hits  = self.hits[(k['id_a'], k['id_b'])];
    T     = self.hit_dendrograms[(k['id_a'], k['id_b'], k['linkage_type'])];

    chrs_a = self.org_chrs[k['id_a']];
    chrs_b = self.org_chrs[k['id_b']];

    C = cluster.calc_clusters(T, hits, chrs_a, chrs_b, alpha=alpha);

    self.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'])] = C;

    return k;
  #edef

  #############################################################################

  def hit_cluster_height(self, k, H):

    if ((k['id_a'], k['id_b']) not in self.hits) or \
       ((k['id_a'], k['id_b']) not in self.hit_distances) or \
       ((k['id_a'], k['id_b'], k['linkage_type']) not in self.hit_dendrograms):
      print "You must run hit_index(), hit_distance() and hit_dendrogram() first!";
      return None;
    #fi

    k['H'] = H;

    hits  = self.hits[(k['id_a'], k['id_b'])];
    T     = self.hit_dendrograms[(k['id_a'], k['id_b'], k['linkage_type'])];
    
    chrs_a = self.org_chrs[k['id_a']];
    chrs_b = self.org_chrs[k['id_b']];
    
    C = cluster.calc_clusters_height(T, hits, chrs_a, chrs_b, H);
    
    self.hit_clusters_height[(k['id_a'], k['id_b'], k['linkage_type'], k['H'])] = C;
    
    return k;
  #edef

  #############################################################################

  def hit_cluster_overlap(self, k):
    """ Identify clusters which are overlapped"""

    if ((k['id_a'], k['id_b'], k['linkage_type'], k['alpha']) not in self.hit_clusters):
      print "You must run hit_cluster() first!";
      return None;
    #fi

    C = self.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'])];

    overlaps = [];
    R        = [];

    for i in xrange(len(C)-1):
      ci  = C[i];
      cio = [];
      for j in xrange(i+1, len(C)):
        cj = C[j];
        rlen = float(ci[2] - ci[1] + 1) / float(cj[2] - ci[1] + 1);
        if ((ci[0] == cj[0]) and
            util.overlap( (ci[1], ci[2]), (cj[1], cj[2])) > 0 and
            (rlen > 0.8 or rlen < 1.25)):
          cio.append(j);
        #fi
      #efor
      if len(cio) > 0:
        overlaps.append(cio + [i]);
      #fi
    #efor

    for (j, cio) in enumerate(overlaps):
      areg_chr   = C[cio[0]][0];
      areg_start = min([ C[i][1] for i in cio ]);
      areg_end   = max([ C[i][2] for i in cio ]);
      areg       = ( self.org_names[k['id_a']], areg_chr, areg_start, areg_end );

      bregs      = [ ( self.org_names[k['id_b']], C[i][3], C[i][4], C[i][5] ) for i in cio ];
      regs       = [ areg ] + bregs;
      R.append( ( str(j), regs));
    #efor

    return R;

   #efor

  #############################################################################

  def hit_clusters_containing(self, k, a_prots, b_prots):


    C  = self.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'])];
    CA = [];
    CB = [];

    for ap in a_prots:
      CA.append([ i for (i,c) in enumerate(C) if ap in c[9] ]);
    #efor
    for bp in b_prots:
      CB.append([ i for (i,c) in enumerate(C) if bp in c[10] ]);
    #efor

    return (CA, CB);
  #edef

  #############################################################################

  def hit_clust_2_reg(self, k, ci, name=None):
    C  = self.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'])][ci];

    #('bri1',   [ ('schco2', 'scaffold_1',  5346259, 5348750), ('agabi', 'scaffold_1',  1865954, 1869064) ] ),

    reg = (name if not(name == None) else str(ci),
           [ ( self.org_names[k['id_a']], C[0], C[1], C[2] ),
             ( self.org_names[k['id_b']], C[3], C[4], C[5] ) ] );

    return reg;
  #edef

  #############################################################################

#eclass

###############################################################################

