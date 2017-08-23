import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
__PC_DIR__ = "%s/pipeline_components" % __INSTALL_DIR__

###############################################################################

import json
dconfig = json.load(open("%s/defaults.json"% __PC_DIR__, "r"))
dconfig.update(config)

###############################################################################

__RUN_DIR__ = os.path.abspath(dconfig["outdir"]) + "/run"
__PROTS_OUTDIR__ = "%s/prots" % __RUN_DIR__
__EXONS_OUTDIR__ = "%s/exons" % __RUN_DIR__
__BLAST_OUTDIR__ = "%s/blast" % __RUN_DIR__
__CLUST_OUTDIR__ = "%s/clust" % __RUN_DIR__
__PROTENY_OUTDIR__ = "%s/proteny" % __RUN_DIR__
__CIRCOS_OUTDIR__ = "%s/circos" % __RUN_DIR__

###############################################################################

rule genProteins:
  input:
    genome = lambda wildcards: dconfig["genomes"][wildcards.genome]["genome"],
    gff    = lambda wildcards: dconfig["genomes"][wildcards.genome]["gff"]
  output:
    fa = "%s/prots.{genome}.fasta"% __PROTS_OUTDIR__
  conda: "%s/conda.yaml"% __PC_DIR__
  shell: """
    gffread -y "{output.fa}.orig" -g "{input.genome}" "{input.gff}"
    sed -e 's/^>\([^ ]\+\).*/>\\1/' {output.fa}.orig \
     | tr '\n>' '\t\n' \
     | sed -e 's/^\([^\t]\+\)\t/>\\1\\n/' \
     | sed -e 's/[.]\([\t]\?\)$/\\1/' \
     | grep -B1 --no-group-separator '^[^>.][^.-]\+$' \
     | tr '\t' '\n' \
     | sed -e '/^$/d' \
     > {output.fa}
  """

###############################################################################

rule exonAminoAcidSequences:
  input:
    prots = lambda wildcards: "%s/prots.%s.fasta" % (__PROTS_OUTDIR__, wildcards.genome),
    gff   = lambda wildcards: dconfig["genomes"][wildcards.genome]["gff"]
  output:
    exonFasta = "%s/exons.{genome}.fasta" % __EXONS_OUTDIR__,
    exonInfo  = "%s/exons.{genome}.info" % __EXONS_OUTDIR__
  run:
    import pipeline_components.gff3utils as Gutils
    import pipeline_components.utils as utils

    G = Gutils.readGFF3File(input.gff)
    F = utils.loadFasta(input.prots)

    exons = utils.indexListBy(G.entries, lambda x: x.attr.get("Parent","__TOPLEVEL__"))

    exonSeqs = {}
    exonInfo = []
    i = 0
    for protein in F:
      if protein not in exons:
        println("'%s' not in GFF file." % protein)
        continue
      #fi

      proteinSeq     = F[protein]
      proteinExons   = sorted([ e for e in exons[protein] if e.type == 'CDS' ], key=lambda e: e.start)
      positiveStrand = proteinExons[0].strand == '+'

      if positiveStrand:
        firstBase = min([ e.start for e in proteinExons])
        lastPosition = 0
        for e in proteinExons:
          exonID     = "G.%s.E.%00000000d" % (wildcards.genome, i)
          i += 1
          exonLength = int((e.end - e.start)/3.0)
          seq        = proteinSeq[lastPosition:lastPosition+exonLength]
          lastPosition += exonLength
          exonSeqs[exonID] = seq
          exonInfo.append((exonID, protein, wildcards.genome, e.seqid, e.start, e.end, e.strand))
        #efor
      else:
        proteinSeq = proteinSeq[::-1]
        firstBase = max([ e.end for e in proteinExons])
        lastPosition = 0
        for e in proteinExons[::-1]:
          exonID     = "G.%s.E.%00000000d" % (wildcards.genome, i)
          i += 1
          exonLength = int((e.end - e.start)/3.0)
          seq        = proteinSeq[lastPosition:lastPosition+exonLength]
          lastPosition += exonLength
          exonSeqs[exonID] = seq
          exonInfo.append((exonID, protein, wildcards.genome, e.seqid, e.start, e.end, e.strand))
        #efor

      #fi

        # Output the FASTA file
      utils.writeFasta(exonSeqs.items(), output.exonFasta)

        # Output the INFO file
      with open(output.exonInfo, "w") as ofd:
        ofd.write("#exonID\tproteinID\tgenome\tseqID\tstart\tend\tstrand\n")
        for ei in exonInfo:
          ofd.write("%s\t%s\t%s\t%s\t%d\t%d\t%s\n" % ei)
        #efor
      #ewith
    

rule exonInfo:
  input:
    exonInfos = expand("%s/exons.{genome}.info" % __EXONS_OUTDIR__, genome=dconfig["genomes"].keys())
  output:
    exonInfo = "%s/exons.info" % __EXONS_OUTDIR__
  shell: """
    cat {input.exonInfos} > {output.exonInfo}
  """

###############################################################################
 # BLAST

rule blastDB:
  input:
    exonFasta = lambda wildcards: "%s/exons.%s.fasta"% (__EXONS_OUTDIR__, wildcards.genome)
  output:
    db = "%s/blastdb.{genome}.db" % __BLAST_OUTDIR__
  conda: "%s/conda.yaml"% __PC_DIR__
  shell: """
    makeblastdb -dbtype prot -in {input.exonFasta} -out {output.db}
    touch {output.db}
  """

from pipeline_components import blastutils as Butils
rule runBlast:
  input:
    db    = lambda wildcards: "%s/blastdb.%s.db" % (__BLAST_OUTDIR__, wildcards.genomeA),
    query = lambda wildcards: "%s/exons.%s.fasta" % (__EXONS_OUTDIR__, wildcards.genomeB)
  output:
    hits = "%s/hits.{genomeA},{genomeB}.tsv" % __BLAST_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  threads: 5
  params:
    blastfields = Butils.blastfields
  shell: """
    blastp -query {input.query} -db {input.db} -outfmt "6 {params.blastfields}" -out {output.hits} -num_threads {threads}
  """

rule reciprocalBlast:
  input:
    hAB = "%s/hits.genome_1,genome_2.tsv" % __BLAST_OUTDIR__,
    hBA = "%s/hits.genome_2,genome_1.tsv" % __BLAST_OUTDIR__
  output:
    hits = "%s/hits.reciprocal.tsv" % __BLAST_OUTDIR__
  run:
    import pipeline_components.blastutils as Butils
    import pipeline_components.utils as utils

    hAB = utils.indexListBy(Butils.readBlastFile(input.hAB), lambda x: (x.sseqid, x.qseqid))
    hBA = utils.indexListBy(Butils.readBlastFile(input.hBA), lambda x: (x.sseqid, x.qseqid))

    recip = []

    for (t1,t2) in hAB.keys():
      if (t2,t1) not in hBA:
        continue
      #fi
      hitsAB = hAB[(t1,t2)]
      hitsBA = hBA[(t2,t1)]

      for h_i in hitsBA:
        for h_j in hitsAB:
          if Butils.blast_hits_overlap(h_i, h_j, lambda x,y: x and y):
            recip.append(h_i)
            break
          #fi
        #efor
      #efor
    #efor

    Butils.writeBlastFile(recip, output.hits)

###############################################################################

 # Recalculating the hit score
rule hitScoreK:
  input:
    hits = rules.reciprocalBlast.output.hits
  output:
    hits = "%s/hits.hitscores.tsv" % __BLAST_OUTDIR__
  run:
    import pipeline_components.blastutils as Butils
    H = Butils.readBlastFile(input.hits)
    augmentedBlastType = Butils.BlastHitType(Butils.augmentedblastfields)

    def hitScoreK(h):
      factor = float((h.send - h.sstart) + (h.qend - h.qstart)) / float(h.slen + h.qlen)
      K = (1.0 - min(1, h.evalue)) * factor
      return K
    #edef

    HK = [ Butils.parseBlastHit(Butils.blastHit2Row(h) + [hitScoreK(h)], augmentedBlastType) for h in H]
    Butils.writeBlastFile(HK, output.hits)

###############################################################################

  # Translate the hits into the genome space
rule genomeSpaceHits:
  input:
    hits = rules.hitScoreK.output.hits,
    exonInfo = rules.exonInfo.output.exonInfo
  output:
    hits = "%s/hits.genomic.tsv"% __BLAST_OUTDIR__
  run:
    import pipeline_components.blastutils as Butils
    import pipeline_components.utils as utils

    hits     = Butils.readBlastFile(input.hits, Butils.augmentedblastfields)
    exonInfo = utils.readExonInfo(input.exonInfo)

    def translateChrPosition(seqid, start, end):
      strand = exonInfo[seqid].strand
      length = exonInfo[seqid].end - exonInfo[seqid].start
      if strand == '+':
        return (start*3 + exonInfo[seqid].start, end*3 + exonInfo[seqid].start)
      else:
        return (exonInfo[seqid].end - (length - end)*3, exonInfo[seqid].end - (length - start)*3)
      #fi
    #edef

    genomicBlastHitType = Butils.BlastHitType(Butils.genomicblastfields)

    ghits = []
    for hit in hits:
      genome1Chr = exonInfo[hit.qseqid].seqid
      genome2Chr = exonInfo[hit.sseqid].seqid

      genome1Region = translateChrPosition(hit.qseqid, hit.qstart, hit.qend)
      genome2Region = translateChrPosition(hit.sseqid, hit.sstart, hit.send)

      ghits.append(Butils.parseBlastHit(hit + (genome1Chr, genome1Region[0], genome1Region[1], genome2Chr, genome2Region[0], genome2Region[1],), genomicBlastHitType ))
    #efor

    Butils.writeBlastFile(ghits, output.hits)

###############################################################################
  # Cluster the hits

def generateChromosomePairs():
  import pipeline_components.gff3utils as Gutils
  GA = Gutils.readGFF3File(dconfig["genomes"]["genome_1"]["gff"])
  GB = Gutils.readGFF3File(dconfig["genomes"]["genome_2"]["gff"])

  pairs = []
  for chrA in GA.seqids:
    for chrB in GB.seqids:
      pairs.append((chrA, chrB))
    #efor
  #efor
  return pairs
#edef

chromosomePairs = generateChromosomePairs()

rule splitHits:
  input:
    hits = rules.genomeSpaceHits.output.hits
  output:
    splitHits = ["%s/chrPairHits.%s.%s.tsv" % (__CLUST_OUTDIR__, chrA, chrB) for (chrA, chrB) in chromosomePairs ]
  run:
    import pipeline_components.blastutils as Butils

    H = Butils.readBlastFile(input.hits, Butils.genomicblastfields)

    chrH = { pair: [] for pair in chromosomePairs }
    for h in H:
      chrH[(h.g1chr,h.g2chr)].append(h)
    #efor

    for (chrA, chrB) in chromosomePairs:
      Butils.writeBlastFile(chrH[(chrA, chrB)], "%s/chrPairHits.%s.%s.tsv" % (__CLUST_OUTDIR__, chrA, chrB))
    #efor

rule makeDendrogram:
  input:
    hits     = lambda wildcards: "%s/chrPairHits.%s.tsv" % (__CLUST_OUTDIR__, wildcards.chrpair),
    exonInfo = rules.exonInfo.output.exonInfo
  output:
    dendrogram = "%s/dendrogram.{chrpair}.tsv" % __CLUST_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  params:
    pc_dir = __PC_DIR__,
  shell: """
    {params.pc_dir}/clusterHits.py "{input.hits}" "{input.exonInfo}" "{output.dendrogram}"
  """

rule dendrograms:
  input:
    dendrograms = expand("%s/dendrogram.{chrpair}.tsv" % __CLUST_OUTDIR__, chrpair = [ '%s.%s' % (a,b) for (a,b) in chromosomePairs ]),
    hits        = expand("%s/chrPairHits.{chrpair}.tsv" % __CLUST_OUTDIR__, chrpair = [ '%s.%s' % (a,b) for (a,b) in chromosomePairs ])
  output:
    fileList = "%s/dataList.tsv" % __CLUST_OUTDIR__
  run:
    import os
    with open(output.fileList, "w") as ofd:
      ofd.write("#chrA\tchrB\thitsfile\tdendrogramfile\n")
      for (chrPair, hitFile, dendrogramFile) in zip(chromosomePairs, input.hits, input.dendrograms):
        if os.stat(hitFile).st_size > 0:
          ofd.write("%s\t%s\t%s\t%s\n" % (chrPair[0], chrPair[1], hitFile, dendrogramFile))
        #fi
      #efor
    #ewith

###############################################################################

rule proteny:
  input:
    dendrogramHitsList = rules.dendrograms.output.fileList,
    exonInfo           = rules.exonInfo.output.exonInfo
  output:
    sigClusters = "%s/proteny.sigClusts.tsv" % __PROTENY_OUTDIR__
  threads: 10
  conda: "%s/conda.yaml"% __PC_DIR__
  params:
    pc_dir = __PC_DIR__,
    pvalue = dconfig["p"],
    cons_thresh = dconfig["cons_thresh"]
  shell: """
    {params.pc_dir}/proteny.py "{input.dendrogramHitsList}" "{input.exonInfo}" "{params.pvalue}" "{params.cons_thresh}" "{threads}" "{output.sigClusters}"
  """

def readSigClusters(sigClustFile):
  from pipeline_components import utils
  C = utils.readColumnFile(sigClustFile, "chrA chrB clusterID startA endA startB endB score pvalue genesA genesB", types="str str str int int int int float float str str")
  return C
#edef

###############################################################################

  # Write the karyotype files
  # All chromosomes are renamed to: "GENOMENAME_SEQID"
rule circosKaryotypeSingle:
  input:
    fasta = lambda wildcards: dconfig["genomes"][wildcards.genome]["genome"]
  output:
    karyotypeFile = "%s/karyotype.{genome}.tsv"% __CIRCOS_OUTDIR__
  run:
    from pipeline_components import utils

    F = utils.loadFasta(input.fasta)

    chrNames = sorted(F.keys())
    nChrs    = len(chrNames)

    with open(output.karyotypeFile, "w") as ofd:
      for i, chrName in enumerate(chrNames):
        color = int(360.0 * (float(i+1)/float(nChrs)));
        ofd.write('chr\t-\t%s_%s\t%s_%s\t0\t%d\thsv(%d,1,1)\n' % (dconfig["genomes"][wildcards.genome]["name"], chrName, dconfig["genomes"][wildcards.genome]["name"], chrName, len(F[chrName]), color))
      #efor
   #ewith

rule circosKaryotypes:
  input:
    karyotypes = expand("%s/karyotype.{genome}.tsv"% __CIRCOS_OUTDIR__, genome=dconfig["genomes"].keys())

rule circosGenes:
  input:
    exonInfo = rules.exonInfo.output.exonInfo
  output:
    genes = "%s/genes.tsv" % __CIRCOS_OUTDIR__
  run:
    from pipeline_components import utils
    E = utils.readExonInfo(input.exonInfo)
    EG = utils.indexListBy(E.values(), lambda e: (e.genome, e.seqid, e.proteinid))
    with open(output.genes, "w") as ofd:
      for (genome, seqid, proteinid) in EG.keys():
        exons = EG[(genome, seqid, proteinid)]
        start = min([ e.start for e in exons ])
        end   = max([ e.end for e in exons ])
        strand = exons[0].strand
        ofd.write('%s_%s\t%d\t%d\tstrand=%s,geneid=%s\n' % (dconfig["genomes"][genome]["name"], seqid, start, end, ('p' if (strand == '+') else 'n'), proteinid, ))
      #efor
    #ewith

rule circosExons:
  input:
    exonInfo = rules.exonInfo.output.exonInfo
  output:
    exons = "%s/exons.tsv" % __CIRCOS_OUTDIR__
  run:
    from pipeline_components import utils
    E = utils.readExonInfo(input.exonInfo)

    with open(output.exons, "w") as ofd:
      for exon in E.values():
        ofd.write('%s_%s\t%d\t%d\tstrand=%s,geneid=%s,exonid=%s\n' % (dconfig["genomes"][exon.genome]["name"], exon.seqid, exon.start, exon.end, ('p' if (exon.strand == '+') else 'n'), exon.proteinid, exon.exonid))
      #efor
    #ewith

rule circosLinks:
  input:
    sigClusters = rules.proteny.output.sigClusters,
  output:
    links = "%s/links.tsv"% __CIRCOS_OUTDIR__
  run:
    import colorsys
    C = readSigClusters(input.sigClusters)

    with open(output.links, "w") as ofd:
      for i, c in enumerate(C):
        #r_int, g_int, b_int = [ int(255*intensity) for intensity in colorsys.hsv_to_rgb(float(color(c))/360.0,1.0,1.0) ]
        r_int, g_int, b_int = [ int(255*intensity) for intensity in colorsys.hsv_to_rgb(1.0,1.0,1.0) ]
        ofd.write('c_%010d\t%s_%s\t%d\t%d\tr_int=%d,g_int=%d,b_int=%d\n' % (i, dconfig["genomes"]["genome_1"]["name"], c.chrA, c.startA, c.endA, r_int, g_int, b_int));
        ofd.write('c_%010d\t%s_%s\t%d\t%d\tr_int=%d,g_int=%d,b_int=%d\n' % (i, dconfig["genomes"]["genome_2"]["name"], c.chrB, c.startB, c.endB, r_int, g_int, b_int));
      #efor
    #ewith

#rule circosBlast:
#  input:
#    hits = rules.genomeSpaceHits.output.hits
#  output:
#    blast = "%s/blast.tsv"% __CIRCOS_OUTDIR__
#  run:
#    from pipeline_components import blastutils as Butils
#    B = Butils.readBlastFile(input.hits, Butils.genomicblastfields)
#
#    with open(output.links, "w") as ofd:
#      for i, c in enumerate(C):
#        #r_int, g_int, b_int = [ int(255*intensity) for intensity in colorsys.hsv_to_rgb(float(color(c))/360.0,1.0,1.0) ]
#        r_int, g_int, b_int = [ int(255*intensity) for intensity in colorsys.hsv_to_rgb(1.0,1.0,1.0) ]
#        odf.write('b_%010d\t%s_%s\t%d\t%d\ta_strand=%s,coverage=%f,score=%f,pident=%f,evalue=%f,bitscore=%f,a_geneid=%s,a_transcriptid=%s,a_exonid=%d,h=%d,s=%0.3f,v=1.0\n'
#        ofd.write('c_%010d\t%s_%s\t%d\t%d\tr_int=%d,g_int=%d,b_int=%d\n' % (i, dconfig["genomes"]["genome_1"]["name"], c.chrA, c.startA, c.endA, r_int, g_int, b_int));
#        ofd.write('c_%010d\t%s_%s\t%d\t%d\tr_int=%d,g_int=%d,b_int=%d\n' % (i, dconfig["genomes"]["genome_2"]["name"], c.chrB, c.startB, c.endB, r_int, g_int, b_int));
#      #efor
#    #ewith

rule circosConf:
  input:
    karyotypes = rules.circosKaryotypes.input.karyotypes,
    genes      = rules.circosGenes.output.genes,
    links      = rules.circosLinks.output.links
  output:
    conf = "%s/configuration.{chr}.conf" % __CIRCOS_OUTDIR__
  run:
    import pipeline_components.circos_core as ccore

    genome2Chromosomes = sorted(list(set([ "%s_%s"% (dconfig["genomes"]["genome_2"]["name"], pair[1]) for pair in chromosomePairs ])))
    chromosomes = "%s_%s;%s" % (dconfig["genomes"]["genome_1"]["name"], wildcards.chr, ';'.join(genome2Chromosomes))

    with open(output.conf, "w") as ofd:
      ofd.write(ccore.makeConf(chromosomes, input.karyotypes, input.genes, input.links))
    #ewith
    

genome1Chromosomes = sorted(list(set([ pair[0] for pair in chromosomePairs ])))
rule circosConfs:
  input:
    confs = expand("%s/configuration.{chr}.conf"% __CIRCOS_OUTDIR__, chr=genome1Chromosomes)
  

rule circosGen:
  input:
    conf = lambda wildcards: "%s/configuration.%s.conf" % (__CIRCOS_OUTDIR__, wildcards.chr)
  output:
    png = "%s/circos.{chr}.png" % __CIRCOS_OUTDIR__,
    svg = "%s/circos.{chr}.svg" % __CIRCOS_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    dirname="`dirname {output.png}`/{wildcards.chr}"
    mkdir -p "$dirname"
    circos -conf {input.conf} -outputdir "$dirname"
    mv "$dirname/circos.png" "{output.png}"
    mv "$dirname/circos.svg" "{output.svg}"
  """

rule circos:
  input:
    figures = expand("%s/circos.{chr}.png" % __CIRCOS_OUTDIR__, chr=genome1Chromosomes)
#rule visualize_data:
#  input:
#    dendrogramHitsList = rules.dendrograms.output.fileList
#    exonInfo           = rules.exonInfo.output.exonInfo,
#    sigClusters        = rules.proteny.output.sigClusters
#  output:
    
