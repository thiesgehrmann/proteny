from ibidas import *;

class circos_data:

  #############################################################################

  def write_data(self, PR, dir):
    self.write_karyotype(PR, dir);
    self.write_genes(PR, dir);
    self.write_exons(PR, dir);
    self.write_blast(PR, dir);
  #edef

  #############################################################################

  def write_blast(self, PR, dir):
    D = dict([(slicename, i) for (i, slicename) in enumerate(PR.__blast_slice_names__)]);
    for (i,j) in PR.blast_hits.keys():
      F = PR.blast_hits[(i,j)];
      filename = '%s/blast.%s.%s.tsv' % (dir, PR.org_names[i], PR.org_names[j])
      fd = open(filename, 'w');

      for (hitid, hit) in enumerate(zip(*F())):
        fd.write('blasthit_%s_%s_%010d\t%s_%s\t%d\t%d\ta_strand=%s,pident=%f,evalue=%f,bitscore=%f,a_geneid=%s,a_transcriptid=%s,a_exonid=%d\n' \
          % (PR.org_names[i], PR.org_names[j], hitid, PR.org_names[i], hit[D['a_chrid']], \
             hit[D['qstart']], hit[D['qend']], \
             ('p' if (hit[D['a_strand']]) == '+' else 'n'), \
             hit[D['pident']], hit[D['evalue']], hit[D['bitscore']], \
             str(hit[D['a_geneid']]), str(hit[D['a_transcriptid']]), hit[D['a_exonid']]));
        fd.write('blasthit_%s_%s_%010d\t%s_%s\t%d\t%d\tb_strand=%s,pident=%f,evalue=%f,bitscore=%f,b_geneid=%s,b_transcriptid=%s,b_exonid=%d\n' \
          % (PR.org_names[i], PR.org_names[j], hitid, PR.org_names[j], hit[D['b_chrid']], \
             hit[D['sstart']], hit[D['send']], \
             ('p' if (hit[D['b_strand']]) == '+' else 'n'), \
             hit[D['pident']], hit[D['evalue']], hit[D['bitscore']], \
             str(hit[D['b_geneid']]), str(hit[D['b_transcriptid']]), hit[D['b_exonid']]));
      #efor
      fd.close();
    #efor
  #edef

  #############################################################################
    # Generate karyotype
  def write_karyotype(self, PR, dir):
    print "Generating karyotype file"
    for (i, name) in enumerate(PR.org_names):
      
      filename = '%s/karyotype.%s.tsv' % (dir, name);
      fd = open(filename, 'w');

      CHRS = zip(*PR.org_genomes[i].Get(_.chrid, _.sequence.Each(lambda x: len(x)))());
      n    = len(CHRS);

      for (j, (id, length)) in enumerate(CHRS):
        fd.write('chr\t-\t%s_%s\t%s_%s\t0\t%d\t%s\n' % (PR.org_names[i], id, PR.org_names[i], id, length, ('hsv(%d,1,1)' % 360.0 * (j / n) )));
      #efor
      fd.close();
    #efor
  #edef

  #############################################################################
    # Generate gene & exon files
  def write_genes(self, PR, dir):

    for (i, E) in enumerate(PR.org_exons):
      filename     = '%s/genes.%s.tab'      % (dir, PR.org_names[i]);
      textfilename = '%s/genes_text.%s.tab' % (dir, PR.org_names[i]);
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
  #edef

  #############################################################################

  def write_exons(self, PR, dir):

    for (i,E) in enumerate(PR.org_exons):
      filename = '%s/exons.%s.tab' % (dir, PR.org_names[i]);
      fd = open(filename, 'w');

      for ex in zip(*E.Get(_.chrid, _.start, _.end, _.strand, _.geneid, _.exonid)()):
        fd.write('%s_%s\t%d\t%d\tstrand=%s,geneid=%s,exonid=%d\n' % (PR.org_names[i], ex[0], ex[1], ex[2], ('p' if (ex[3] == '+') else 'n'), str(ex[4]), ex[5]));
      #efor
      fd.close();
    #efor
  #edef

  #############################################################################

#eclass

