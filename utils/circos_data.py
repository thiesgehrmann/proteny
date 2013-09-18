from ibidas import *;

#############################################################################

def write_data(PR, dir):
  files = {};

  files['karyotype'], chr_colors = write_karyotype(PR, dir);
  files['genes']                 = write_genes(PR, dir);
  files['exons']                 = write_exons(PR, dir);
  files['blast']                 = write_blast(PR, dir, chr_colors);
  files['clusts']                = write_clusts(PR, dir, chr_colors);

  return files;
#edef

#############################################################################

def write_clusts(PR, dir, chr_colors):
  fnames = [];
  for (i, j) in PR.hit_clusters.keys():
    (HD, CS) = PR.hit_clusters[(i,j)];

    color_org = (PR.org_names[i], 0) if len(chr_colors[PR.org_names[i]]) < len(chr_colors[PR.org_names[j]]) else (PR.org_names[j], 4);
    color = lambda h: chr_colors[color_org[0]][h[color_org[1]]];

    for (l, cs) in enumerate(CS):
      hd = HD[l];
      filename = '%s/clusters_%s_%s_%d_avg_sum_subset.tsv' % (dir, PR.org_names[i], PR.org_names[j], cs);
      fnames.append(filename);
      fd = open(filename, 'w');
      for (k, h) in enumerate(hd):
        fd.write('clusters_%s_%s_%010d\t%s_%s\t%d\t%d\ta_chrid=%d,nhits=%d,score=%f,prots_a=%s,h=%d\n' % (PR.org_names[i], PR.org_names[j], k, PR.org_names[i], h[0], h[2], h[3], h[1], h[8], h[9], ';'.join(['%s' % id for id in h[10]]), color(h)));
        fd.write('clusters_%s_%s_%010d\t%s_%s\t%d\t%d\tb_chrid=%d,nhits=%d,score=%f,prots_b=%s,h=%d\n' % (PR.org_names[i], PR.org_names[j], k, PR.org_names[j], h[4], h[6], h[7], h[5], h[8], h[9], ';'.join(['%s' % id for id in h[11]]), color(h)));
      #efor
      fd.close();
    #efor
  #efor
  return fnames;
#edef

#############################################################################

def write_blast(PR, dir, chr_colors):
  fnames = [];
  D = dict([(slicename, i) for (i, slicename) in enumerate(PR.__blast_slice_names__)]);
  for (i,j) in PR.blast_hits.keys():
    F = PR.blast_hits[(i,j)];
    filename = '%s/blast.%s.%s.tsv' % (dir, PR.org_names[i], PR.org_names[j]);
    fnames.append(filename);
    fd = open(filename, 'w');

    color_org = (PR.org_names[i], 'a_chrid') if len(chr_colors[PR.org_names[i]]) < len(chr_colors[PR.org_names[j]]) else (PR.org_names[j], 'b_chrid');
    color = lambda h: chr_colors[color_org[0]][h[D[color_org[1]]]];

    H = zip(*F());
    for (hitid, hit) in enumerate(H):
      fd.write('blasthit_%s_%s_%010d\t%s_%s\t%d\t%d\ta_strand=%s,pident=%f,evalue=%f,bitscore=%f,a_geneid=%s,a_transcriptid=%s,a_exonid=%d,h=%d,s=%0.3f,v=1.0\n' \
        % (PR.org_names[i], PR.org_names[j], hitid, PR.org_names[i], hit[D['a_chrid']], \
           hit[D['a_start']], hit[D['a_end']], \
           ('p' if (hit[D['a_strand']]) == '+' else 'n'), \
           hit[D['pident']], hit[D['evalue']], hit[D['bitscore']], \
           str(hit[D['a_geneid']]), str(hit[D['a_transcriptid']]), hit[D['a_exonid']], \
           color(hit), hit[D['bitscore']]));
      fd.write('blasthit_%s_%s_%010d\t%s_%s\t%d\t%d\tb_strand=%s,pident=%f,evalue=%f,bitscore=%f,b_geneid=%s,b_transcriptid=%s,b_exonid=%d,h=%d,s=%0.3f,v=1.0\n' \
        % (PR.org_names[i], PR.org_names[j], hitid, PR.org_names[j], hit[D['b_chrid']], \
           hit[D['b_start']], hit[D['b_end']], \
           ('p' if (hit[D['b_strand']]) == '+' else 'n'), \
           hit[D['pident']], hit[D['evalue']], hit[D['bitscore']], \
           str(hit[D['b_geneid']]), str(hit[D['b_transcriptid']]), hit[D['b_exonid']], \
           color(hit), hit[D['bitscore']]));
    #efor
    fd.close();
  #efor
  return fnames;
#edef

#############################################################################
  # Generate karyotype
def write_karyotype(PR, dir):
  fnames = [];
  chr_colors = {};

  f = lambda x, y: y + 1 -(x + (1 - (x % 2))*(y/2 - x/2) - (x%2)*(x/2));
  f = lambda x, y: ((((x-1)*(y/2 -1)) % y) + 1);
  f = lambda x, y: x;

  for (i, name) in enumerate(PR.org_names):
    
    filename = '%s/karyotype.%s.tsv' % (dir, name);
    fnames.append(filename);
    fd = open(filename, 'w');

    CHRS = zip(*PR.org_genomes[i].Get(_.chrid, _.sequence.Each(lambda x: len(x)))());
    n    = len(CHRS);
    #n = n - (n%2);

    for (j, (id, length)) in enumerate(CHRS):
      color = int(360.0 * (float(f(j+1,n))/float(n)));
      fd.write('chr\t-\t%s_%s\t%s_%s\t0\t%d\thsv(%d,1,1)\n' % (PR.org_names[i], id, PR.org_names[i], id, length, color));
      if PR.org_names[i] not in chr_colors:
        chr_colors[PR.org_names[i]] = {};
      #fi
      chr_colors[PR.org_names[i]][id] = color;
    #efor
    fd.close();
  #efor

  return (fnames, chr_colors);
#edef

#############################################################################
  # Generate gene & exon files
def write_genes(PR, dir):
  fnames = { 'text':[], 'data':[] };


  for (i, E) in enumerate(PR.org_exons):
    filename     = '%s/genes.%s.tsv'      % (dir, PR.org_names[i]);
    textfilename = '%s/genes_text.%s.tsv' % (dir, PR.org_names[i]);

    fnames['data'].append(filename);
    fnames['text'].append(textfilename);

    fd     = open(filename, 'w');
    textfd = open(textfilename, 'w');

    G = E.Without(_.sequence).GroupBy(_.geneid).Get(_.chrid[0], Min(_.start), Max(_.end), _.strand[0], _.geneid);

    for gene in zip(*G()):
      fd.write(    '%s_%s\t%d\t%d\tstrand=%s,geneid=%s\n'     % (PR.org_names[i], gene[0], gene[1], gene[2], ('p' if (gene[3] == '+') else 'n'), str(gene[4]), ));
      textfd.write('%s_%s\t%d\t%d\t%s\tstrand=%s,geneid=%s\n' % (PR.org_names[i], gene[0], gene[1], gene[2], str(gene[4]), ('p' if (gene[3] == '+') else 'n'), str(gene[4])));
    #efor
    fd.close();
    textfd.close();
  #efor

  return fnames;
#edef

#############################################################################

def write_exons(PR, dir):
  fnames = [];
  for (i,E) in enumerate(PR.org_exons):
    filename = '%s/exons.%s.tsv' % (dir, PR.org_names[i]);
    fnames.append(filename);
    fd = open(filename, 'w');

    for ex in zip(*E.Get(_.chrid, _.start, _.end, _.strand, _.geneid, _.exonid)()):
      fd.write('%s_%s\t%d\t%d\tstrand=%s,geneid=%s,exonid=%d\n' % (PR.org_names[i], ex[0], ex[1], ex[2], ('p' if (ex[3] == '+') else 'n'), str(ex[4]), ex[5]));
    #efor
    fd.close();
  #efor

  return fnames;
#edef

#############################################################################


