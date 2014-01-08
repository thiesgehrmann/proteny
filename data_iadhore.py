#!/usr/bin/python

from ibidas import *;
import sys;
import os, errno;
sys.path.append('./utils');
import seq as sequtils;

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

#F = pd.make_ini(sc, ag, '/home/nfs/thiesgehrmann/groups/w/phd/tasks/synteny/i-adhore-3.0.01/build/experiments/schco2,agabi2', name1="schco2", name2="agabi2");

