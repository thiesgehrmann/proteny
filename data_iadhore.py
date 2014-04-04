#!/usr/bin/python

from ibidas import *;
import sys;
import os, errno;
sys.path.append('./utils');
import seq as sequtils;

from visualization import circos_data as cdata;
reload(cdata);
from visualization import circos_region as cr;
from visualization import circos_chr    as cc;

from utils import cluster;
reload(cluster);
from utils import cluster_null as null;

import bisect;

###############################################################################

td = "/home/nfs/thiesgehrmann/groups/w/phd/tasks/synteny/proteny";
bd = "/tudelft.net/staff-groups/ewi/insy/DBL/thiesgehrmann/w/phd/ibidas_blasts";

###############################################################################

#F = data_iadhore.make_ini(data.aniger_n402(), data.aniger_513_88(), '/home/nfs/thiesgehrmann/groups/w/phd/tasks/synteny/i-adhore-3.0.01/build/experiments/aspni')
#F = data_iadhore.make_ini(data.schco2(), data.agabi2(), '/home/nfs/thiesgehrmann/groups/w/phd/tasks/synteny/i-adhore-3.0.01/build/experiments/schco2,agabi2')

###############################################################################

def prepare_data(org1, org2, dir='./'):

  name1, genes1, genome1 = org1;
  name2, genes2, genome2 = org2;

  org1_g = genes1.Unique(_.transcriptid).Sort(_.chrid, _.start).Detect().Get(_.chrid, _.transcriptid, _.transcriptid.Cast(str) + _.strand);
  org2_g = genes2.Unique(_.transcriptid).Sort(_.chrid, _.start).Detect().Get(_.chrid, _.transcriptid, _.transcriptid.Cast(str) + _.strand);

  #org1_g = org1.genes[_.region_type == 'exon'].To(_.attribute, Do=_.Each(_.split(' ')[-1])).Unique(_.attribute).Sort(_.assemblyid, _.start).Detect().Get(_.assemblyid, _.attribute, _.attribute.Cast(str) + _.strand).Copy();
  #org2_g = org2.genes[_.region_type == 'exon'].To(_.attribute, Do=_.Each(_.split(' ')[-1])).Unique(_.attribute).Sort(_.assemblyid, _.start).Detect().Get(_.assemblyid, _.attribute, _.attribute.Cast(str) + _.strand).Copy();

  org1_p = translate_from_sources(genome1, genes1);
  org2_p = translate_from_sources(genome2, genes2);
  #org1_p = org1.prots#.To(_.proteinid, Do=_.TakeFrom(org1_g.Get(_.attribute, _.result)));
  #org2_p = org2.prots#.To(_.proteinid, Do=_.TakeFrom(org2_g.Get(_.attribute, _.result)));

  B = (org1_p | Blast(reciprocal=True, folder=bd) | org2_p).Copy();
  #B = B.Unique(_.L.result, _.R.result);
  #B = B.Get(name1 + '_' +_.L.result, name2 + '_' + _.R.result);
  B = B.Unique(_.L.transcriptid, _.R.transcriptid);
  B = B.To(_.L.transcriptid, Do=_.Cast(str));
  B = B.To(_.R.transcriptid, Do=_.Cast(str));
  B = B.Get(name1 + '_' +_.L.transcriptid, name2 + '_' + _.R.transcriptid);

  mkdir_p(dir);
  mkdir_p('%s/%s_lst' % (dir, name1));
  mkdir_p('%s/%s_lst' % (dir, name2));

  B_loc = '%s_%s_blasttable.tsv' % (name1, name2);
  Export(B, '%s/%s' % (dir, B_loc), names=False);

  # Make directory
  # Export for each chromosome.
  lsts1 = []
  for chr in org1_g.chrid.Unique().Sort()():
    C_loc = '%s_lst/%s.tsv' % (name1, str(chr));
    Export(org1_g[_.chrid == chr].Get(name1 + '_' + _.result), '%s/%s' % (dir, C_loc), names=False);
    lsts1.append(C_loc);
  #efor
  lsts2 = [];
  for chr in org2_g.chrid.Unique().Sort()():
    C_loc = '%s_lst/%s.tsv' % (name2, str(chr));
    Export(org2_g[_.chrid == chr].Get(name2 + '_' + _.result), '%s/%s' % (dir, C_loc), names=False);
    lsts2.append(C_loc);
  #efor

  return (lsts1, lsts2, B_loc);
#edef

###############################################################################

def write_ini(lists1, lists2, blast_table, dir, name1='org1', name2='org2'):

  fd = open('%s/experiment.ini' % dir, 'w');

  fd.write("genome= %s\n" % name1);
  for (i, chr) in enumerate(lists1):
    fd.write("%d %s\n" % (i, chr));
  #efor
  fd.write('\n');
  fd.write("genome= %s\n" % name2);
  for (i, chr) in enumerate(lists2):
    fd.write("%d %s\n" % (i, chr));
  #efor
  fd.write('\n');
  fd.write('blast_table= %s\n' % blast_table);
  fd.write('output_path= results/\n');
  fd.write('write_stats= true\n');
  fd.write('table_type= pairs\n');
  fd.write('\n');
  fd.write('cluster_type= collinear\n');
  fd.write('gap_size= 55\n');
  fd.write('cloud_gap_size= 20\n');
  fd.write('cloud_cluster_gap= 25\n');
  fd.write('cloud_filter_method= binomial\n');
  fd.write('max_gaps_in_alignment= 35\n');
  fd.write('cluster_gap= 55\n');
  fd.write('tandem_gap= 17\n');
  fd.write('\n');
  fd.write('q_value= 0.75\n');
  fd.write('prob_cutoff= 0.01\n');
  fd.write('anchor_points= 3\n');
  fd.write('alignment_method= nw\n');
  fd.write('level_2_only= false\n');
  fd.write('multiple_hypothesis_correction= FDR\n');

  fd.close();

  return '%s/experiment.ini' % dir;
#edef

###############################################################################

def make_ini(org1, org2, dir):
  (g1, g2, B) = prepare_data(org1, org2, dir=dir);
  return write_ini(g1, g2, B, dir, name1=org1[0], name2=org2[0]);
#edef

###############################################################################

def mkdir_p(path):
    # By lzot from http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
#edef

###############################################################################

def translate_from_sources(genome, genes):
  
  EG = genes.GroupBy(_.transcriptid).Get(_.chrid[0], _.start, _.end, _.strand[0], _.geneid[0], _.transcriptid, _.exonid).Copy();
  
  dna = dict(zip(*genome()));
  
  prots = [];
  
  for gene in zip(*EG()):
    aid, starts, ends, strand, geneid, transcriptid, exonids = gene;
    chrdna = dna[aid];
    transcript = '';
    exonsource = [0];
    for (s, e, id) in sorted(zip(starts, ends, exonids), key=lambda x: x[2], reverse=False):
      tr_ex       = chrdna[s-1:e];
      transcript += sequtils.revcomp(tr_ex) if (strand == '-') else tr_ex;
      exonsource += [exonsource[-1] + int((e-(s-1))/3)];
    #efor
    
    aa_seq = sequtils.translate(transcript);
    prots.append((aid, strand, geneid, transcriptid, aa_seq));

  #efor
  return Rep(prots) / ('chrid', 'strand', 'geneid', 'transcriptid', 'sequence');

#edef

###############################################################################

class iadhore_fit_into_PR:

  org_names    = [];
  org_genomes  = [];

  hit_clusters = {};

  k = {
  'id_a'         : 0,
  'id_b'         : 1,
  'linkage_type' : 'NA',
  'alpha'        : 0.0,
  'cut'          : 'NA',
  'nd'           : 'NA' }

  def __init__(self, org1, org2, hit_clusters):
    k = self.k;

    self.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd'])] = hit_clusters;
    self.org_names   = [ org1[0], org2[0] ];
    self.org_genomes = [ org1[2], org2[2] ];
    self.org_genes   = [ org1[1], org2[1] ];
    self.org_exons   = [ x.AddSlice('sequence', x.Get(0)) for x in self.org_genes ];
    self.blast_hits  = {};

  #edef

  def key_s(self, k):
    names  = '_'.join([self.org_names[k['id_a']], self.org_names[k['id_b']]]);
    params = '_'.join([k['linkage_type'] ]);
    clust  = '_'.join([k['cut'], k['nd'], str(k['alpha'])]);
    
    return '_'.join([names, params, clust]);
  #edef

#eclass

###############################################################################

def iadhore_results_circos(multiplicon_pairs, org1, org2, order_weird=False):

  genes1 = org1[1];
  genes2 = org2[1];
  genes1 = genes1.Get(_.chrid, _.geneid, _.start, _.end).GroupBy(_.geneid).Get(_.chrid[0], _.geneid, _.start.Min(), _.end.Max()) / ('c1', 'g1', 's1', 'e1');
  genes2 = genes2.Get(_.chrid, _.geneid, _.start, _.end).GroupBy(_.geneid).Get(_.chrid[0], _.geneid, _.start.Min(), _.end.Max()) / ('c2', 'g2', 's2', 'e2');

  D = Read(multiplicon_pairs);
  D = D.Get(_.f1, _.f2, _.f3);

  ND = [];
  for d in zip(*D()):
    clustid = d[0];
    gene_1  = d[1];
    gene_2  = d[2];
    if gene_1.split('_')[0] == org1[0]:
      L = (clustid, gene_1, gene_2);
    else:
      L = (clustid, gene_2, gene_1);
    #fi
    ND.append(L);
  #fi
  D = Rep(ND);
  D = D.To(_.f1, Do=_.Each(lambda x: x.split('_')[1]));
  D = D.To(_.f2, Do=_.Each(lambda x: x.split('_')[1])).Detect();
  D = D / ('clusterid', 'gene_1', 'gene_2');

  Z = (D | Match(_.gene_1, _.g1, merge_same='equi') | genes1) | Match(_.gene_2, _.g2, merge_same='equi') | genes2;
  Z = Z.GroupBy(_.clusterid);
 
  Z = Z.To(_.c1, Do=_[0]);
  Z = Z.To(_.c2, Do=_[0]); 
  Z = Z.To(_.s1, Do=_.Min());
  Z = Z.To(_.s2, Do=_.Min());
  Z = Z.To(_.e1, Do=_.Max());
  Z = Z.To(_.e2, Do=_.Max());

  clusts = [];
  for z in zip(*Z()):
    cd = (z[3], z[4], z[5], z[6], z[7], z[8], 0.0, (1,0,0,0), 0.0, set(z[1]), set(z[2]));
    clusts.append(cd);
  #efor

  PR = iadhore_fit_into_PR(org1, org2, clusts);

  return PR;

#edef

###############################################################################

def score_iadhore_PR(PR, PR_full):

  K = PR.k;

  hits  = PR_full.hits[(K['id_a'], K['id_b'])];

  hit_index     = {};
  hit_index_rev = {};
  for h in hits:
    k = (str(h[1]),str(h[4]));
    if k not in hit_index:
      hit_index[k] = [];
    #fi
    hit_index[k].append(h);
  #efor
  hit_index_rev 
  for k in hit_index.keys():
    hi               = sorted(hit_index[k], key=lambda x: x[7]);
    hir              = sorted(hit_index[k], key=lambda x: x[9])
    hit_index[k]     = ([x[7] for x in hi],  [x[0] for x in hi]);
    hit_index_rev[k] = ([x[9] for x in hir], [x[0] for x in hir]);
  #efor
  
  chrs_a = PR_full.org_chrs[K['id_a']];
  chrs_b = PR_full.org_chrs[K['id_b']];
  
    # Precompute max hit scores per exon
  hitA = {};
  hitB = {};
  v     = [ (2,3,hitA), (5,6,hitB) ];
  
  for (i, h) in enumerate(hits):
    for (j, k, H) in v:
      k = (h[j], h[k]);
      if not(k in H):
        H[k] = i;
      elif h[12] > hits[H[k]][12]:
         H[k] = i;
      #fi
    #efor
  #efor
  
  #T = dict([T.items()[0]]);

  tests = 0;
            # Chromosome 1, chromosome2, index, score, pvalue)
  #start_clusts = [ (k[0], k[1], -1, None, None) for k in T.keys() ];
  clusts = PR.hit_clusters[(K['id_a'], K['id_b'], K['linkage_type'], K['alpha'], K['cut'], K['nd'])]
  
  nd = cluster.null_dist(dist=null.cluster_null_score_strict_smart,
                 tests = 1,
                 scores=[ h[12] for h in hits],
                 hits=hits, hitA=hitA, hitB=hitB,
                 sperm=1000,
                 alpha = 0.05,
                 storage='null_dist_store_cmp_clt.dat');
  nd.update_mperm(len(clusts));

  NC = [];
  for C in clusts:
      # Get the hits within the region defined by the cluster
    k = (str(C[0]), str(C[3]));
    his, hi   = hit_index[k];
    hirs, hir = hit_index_rev[k];

    H_a_s = bisect.bisect_left(his, C[1]);
    H_a_e = bisect.bisect_left(his[H_a_s:], C[2]);
    H_b_s = bisect.bisect_left(hirs, C[4]);
    H_b_e = bisect.bisect_left(hirs[H_b_s:], C[5]);
    H_a   = set(hi[H_a_s:H_a_s+H_a_e]);
    H_b   = set(hir[H_b_s:H_b_s+H_b_e]);
    H_i = list(H_a & H_b);
    H   = [ hits[i] for i in H_i ];

    if len(H) == 0:
      C_info = (0, 0, 0, 0, 0);
      pvalue = 1.0;
    else:
      C_info = cluster.score_hits(H, hits, chrs_a, chrs_b, hitA, hitB);
      C_info = (C_info[0], len(H), len(C_info[5]), len(C_info[6]));
      score, n, nz_nue_a, nz_nue_b = C_info;
    #fi

    (zscore, pvalue) = nd.pvalue(score=score, n=n, nue_a=nz_nue_a, nue_b=nz_nue_b);

    nc = (int(C[0]), C[1], C[2], int(C[3]), C[4], C[5], C[6], C_info, pvalue, C[9], C[10], H_i);
    NC.append(nc);

  #efor

  PR.hit_clusters[(K['id_a'], K['id_b'], K['linkage_type'], K['alpha'], K['cut'], K['nd'])] = NC;

  return NC;
#edef

#edef

###############################################################################

def viz_iadhore_PR(PR, circdir):

  k = PR.k;

  files = cdata.write_data(PR, k, '%s/data' % circdir);
  KS    = PR.key_s(k);

  bchrs = ["%s_%s" % (PR.org_names[k['id_b']], str(i)) for i in PR.org_genomes[k['id_b']].Get(0)() ];
  for i in xrange(PR.org_genomes[k['id_a']].Shape()()):
    cc.circos_chr(files, bchrs + ["%s_%d" % (PR.org_names[k['id_a']], i+1)], ["%s_%d=0.4r" % (PR.org_names[k['id_a']], i+1)], circdir, "scaffold_%02d_%s" % (i+1, KS) );
    print "scaffold_%02d_%s" % (i+1, KS);
  #efor

#edef

###############################################################################

  #(a_chr, a_start, a_end, b_chr, b_start, b_end, n_hits, score, prots_a, prots_b)




#F   = pd.make_ini(sc, ag, '/home/nfs/thiesgehrmann/groups/w/phd/tasks/synteny/i-adhore-3.0.01/build/experiments/schco2,agabi2', name1="schco2", name2="agabi2");

# sc = data.schco2();
# ag = data.agabi2();
# PRi_basid = idata.iadhore_results_circos('~/groups/w/phd/tasks/synteny/i-adhore-3.0.01/build/experiments/schco2,agabi2/results/multiplicon_pairs.txt', sc, ag)
# idata.score_iadhore_PR(PRi_basid, PR_basid)
# idata.viz_iadhore_PR(PRi_basid, 'visualizations/iadhore/basid');

# an1 = data.aniger_n402();
# an2 = data.aniger_513_88();
# PRi_aniger = idata.iadhore_results_circos('~/groups/w/phd/tasks/synteny/i-adhore-3.0.01/build/experiments/aspni/results/multiplicon_pairs.txt', an1, an2)
# idata.score_iadhore_PR(PRi_aniger, PR_aspni);
# idata.viz_iadhore_PR(PRi_aniger, 'visualizations/iadhore/aniger');


#PRi = iadhore_results_circos('multiplicon_pairs.txt', org1, org2);
