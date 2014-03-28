from ibidas import *;
import os;
import errno;
import colorsys;

#############################################################################

def write_data(PR, k, dir):
  """Write data for task k in directory dir"""
  files = {};

  mkdir_p(dir);

  files['karyotype'], chr_colors = write_karyotype(PR, k, dir);
  files['genes']                 = write_genes(PR, k, dir);
  files['exons']                 = write_exons(PR, k, dir);
  files['blast']                 = write_blast(PR, k, dir, chr_colors);
  files['clusts']                = write_clusts(PR, k, dir, chr_colors);

  return files;
#edef

#############################################################################

def write_clusts(PR, k, dir, chr_colors):
  fnames = [];

  i = k['id_a'];
  j = k['id_b'];

  C = PR.hit_clusters[(k['id_a'], k['id_b'], k['linkage_type'], k['alpha'], k['cut'], k['nd'])];

  color_org = (PR.org_names[i], 0) if len(chr_colors[PR.org_names[i]]) < len(chr_colors[PR.org_names[j]]) else (PR.org_names[j], 3);
  color = lambda h: chr_colors[color_org[0]][str(h[color_org[1]])];

  filename = '%s/clusters_%s.tsv' % (dir, PR.key_s(k));
  fnames.append(filename);
  fd = open(filename, 'w');
  for (k, c) in enumerate(C):
    # 0      1        2      3      4        5      6       7      8        9
    #(a_chr, a_start, a_end, b_chr, b_start, b_end, n_hits, score, prots_a, prots_b)
    r_int, g_int, b_int = [ int(255*intensity) for intensity in colorsys.hsv_to_rgb(float(color(c))/360.0,1.0,1.0) ]
    fd.write('clusters_%s_%s_%010d\t%s_%s\t%d\t%d\tnhits=%d,score=%f,n=%d,nue_a=%d,nue_b=%d,p=%f,prots_a=%s,h=%d,r_int=%d,g_int=%d,b_int=%d\n' % (PR.org_names[i], PR.org_names[j], k, PR.org_names[i], c[0], c[1], c[2], c[6], c[7][0],c[7][1],c[7][2],c[7][3], c[8], ';'.join([str(id) for id in c[9]]),  color(c), r_int, g_int, b_int));
    fd.write('clusters_%s_%s_%010d\t%s_%s\t%d\t%d\tnhits=%d,score=%f,n=%d,nue_a=%d,nue_b=%d,p=%f,prots_b=%s,h=%d,r_int=%d,g_int=%d,b_int=%d\n' % (PR.org_names[i], PR.org_names[j], k, PR.org_names[j], c[3], c[4], c[5], c[6], c[7][0],c[7][1],c[7][2],c[7][3], c[8], ';'.join([str(id) for id in c[10]]), color(c), r_int, g_int, b_int));
  #efor
  fd.close();

  return fnames;
#edef

#############################################################################

def write_blast(PR, k, dir, chr_colors):

  fnames = [];
  i      = k['id_a'];
  j      = k['id_b'];

  if (i,j) not in PR.blast_hits:
   return fnames;
  #fi

  D = dict([(slicename, indx) for (indx, slicename) in enumerate(PR.__blast_slice_names__)]);

  F = PR.blast_hits[(i,j)];
  filename = '%s/blast.%s.%s.tsv' % (dir, PR.org_names[i], PR.org_names[j]);
  fnames.append(filename);
  fd = open(filename, 'w');

  color_org = (PR.org_names[i], 'a_chrid') if len(chr_colors[PR.org_names[i]]) < len(chr_colors[PR.org_names[j]]) else (PR.org_names[j], 'b_chrid');
  color = lambda h: chr_colors[color_org[0]][str(h[D[color_org[1]]])];

  H = zip(*F());
  for (hitid, hit) in enumerate(H):
    fd.write('blasthit_%s_%s_%010d\t%s_%s\t%d\t%d\ta_strand=%s,coverage=%f,score=%f,pident=%f,evalue=%f,bitscore=%f,a_geneid=%s,a_transcriptid=%s,a_exonid=%d,h=%d,s=%0.3f,v=1.0\n' \
      % (PR.org_names[i], PR.org_names[j], hitid, PR.org_names[i], hit[D['a_chrid']], \
         hit[D['a_start']], hit[D['a_end']], \
         ('p' if (hit[D['a_strand']]) == '+' else 'n'), \
         hit[D['coverage']], hit[D['score']], hit[D['pident']], hit[D['evalue']], hit[D['bitscore']], \
         str(hit[D['a_geneid']]), str(hit[D['a_transcriptid']]), hit[D['a_exonid']], \
         color(hit), hit[D['bitscore']]));
    fd.write('blasthit_%s_%s_%010d\t%s_%s\t%d\t%d\tb_strand=%s,coverage=%f,score=%f,pident=%f,evalue=%f,bitscore=%f,b_geneid=%s,b_transcriptid=%s,b_exonid=%d,h=%d,s=%0.3f,v=1.0\n' \
      % (PR.org_names[i], PR.org_names[j], hitid, PR.org_names[j], hit[D['b_chrid']], \
         hit[D['b_start']], hit[D['b_end']], \
         ('p' if (hit[D['b_strand']]) == '+' else 'n'), \
         hit[D['coverage']], hit[D['score']], hit[D['pident']], hit[D['evalue']], hit[D['bitscore']], \
         str(hit[D['b_geneid']]), str(hit[D['b_transcriptid']]), hit[D['b_exonid']], \
         color(hit), hit[D['bitscore']]));
  #efor
  fd.close();

  return fnames;
#edef

#############################################################################
  # Generate karyotype
def write_karyotype(PR, k, dir):
  fnames = [];
  chr_colors = {};

  for i in [k['id_a'], k['id_b']]:
    name     = PR.org_names[i];
    filename = '%s/karyotype.%s.tsv' % (dir, name);
    fnames.append(filename);
    fd = open(filename, 'w');

    CHRS = zip(*PR.org_genomes[i].Get(_.chrid, _.sequence.Each(lambda x: len(x)))());
    n    = len(CHRS);

    for (j, (id, length)) in enumerate(CHRS):
      color = int(360.0 * (float(j+1)/float(n)));
      fd.write('chr\t-\t%s_%s\t%s_%s\t0\t%d\thsv(%d,1,1)\n' % (PR.org_names[i], id, PR.org_names[i], id, length, color));
      if PR.org_names[i] not in chr_colors:
        chr_colors[PR.org_names[i]] = {};
      #fi
      chr_colors[PR.org_names[i]][str(id)] = color;
    #efor
    fd.close();
  #efor

  return (fnames, chr_colors);
#edef

#############################################################################
  # Generate gene files
def write_genes(PR, k, dir):
  fnames = { 'text':[], 'data':[] };

  for i in [k['id_a'], k['id_b']]:
    E            = PR.org_exons[i];
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

def write_exons(PR, k, dir):
  fnames = [];
  for i in [k['id_a'], k['id_b']]:
    E        = PR.org_exons[i];
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

def mkdir_p(path):
    # By lzot from http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
#edef

