#!/usr/bin/env nextflow

species = [
  'Salmonella-enterica': [
    'referenceURL': 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz',
    'readsID': '2492428'
  ],
  // 'Staphylococcus-aureus': [
  //   'referenceURL': 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000013425.1_ASM1342v1/GCA_000013425.1_ASM1342v1_genomic.fna.gz',
  //   'readsID': '1274026'
  // ]
]

BASE_PATH = System.getProperty('user.dir') + '/'
BIN = BASE_PATH + '../bin'
ADAPTERS = BASE_PATH + '../adapters'

THREADS = 20

echo true

// === DOWNLOAD ===

process downloadReference {
  container true

  input: val referenceURL from species.collect { it.value.referenceURL }
  output: file 'reference.genomic.fna.gz' into referenceGenomeGz

  """
  appropriate/curl $referenceURL -o reference.genomic.fna.gz
  """
}

// fork reference genome into three other channels
( referenceGenomeGz1,
  referenceGenomeGz2,
  referenceGenomeGz3,
  referenceGenomeGz4 ) = referenceGenomeGz.into(4)

process downloadSRA {
  container 'bionode/bionode-ncbi'

  input: val readsID from species.collect { it.value.readsID }
  output: file '**/*.sra' into reads

  """
  bionode-ncbi download sra $readsID > tmp
  """
}

// === EXTRACT/DECOMPRESS ===

process extractSRA {
  container 'inutano/sra-toolkit'

  input: file read from reads
  output: file '*.fastq.gz' into samples

  """
  fastq-dump --split-files --skip-technical --gzip $read
  """
}

( samples1,
  samples2 ) = samples.into(2)

process decompressReference {
  container 'biodckrdev/htslib'

  input: file referenceGenome from referenceGenomeGz1
  output: file 'reference.genomic.fna' into referenceGenomes

  """
  bgzip -d $referenceGenome --stdout > reference.genomic.fna
  """
}

( referenceGenomes1,
  referenceGenomes2 ) = referenceGenomes.into(2)

// === TRIMMING ===

process trim {
  container 'thejmazz/ployglot-ngs-01'

  input: file reads from samples1
  output: file '*pe.fastq.gz' into trimmomaticDump

  """
  java -jar $BIN/trimmomatic-0.36.jar PE -phred33 \
  $reads \
  reads_1.trim.pe.fastq.gz reads_1.trim.se.fastq.gz \
  reads_2.trim.pe.fastq.gz reads_2.trim.se.fastq.gz \
  ILLUMINACLIP:$ADAPTERS/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
  """
}

process mergeTrimEnds {
  container 'biodckr/bwa'

  input: file ends from trimmomaticDump
  output: file 'reads.trim.pe.fastq.gz' into mergedTrim

  """
  seqtk mergepe $ends | gzip > reads.trim.pe.fastq.gz
  """
}

// === MAPPING ===

process indexReference {
  container 'biodckr/bwa'

  input: file reference from referenceGenomeGz2
  output: file '*.gz.*' into referenceIndexes

  """
  bwa index $reference
  """
}

( referenceIndexes1,
  referenceIndexes2 ) = referenceIndexes.into(2)

// === BWA MEM ===

process alignReadsVanilla {
  container 'bwa-samtools-bcftools'

  input:
    file reference from referenceGenomeGz3
    file referenceIndex from referenceIndexes1
    file sample from samples2
  output: file 'reads.sam' into readsUnsorted


  """
  bwa mem -t $THREADS $reference $sample | samtools view -Sbh - -o reads.sam
  """
}

process alignReads_trim {
  container 'bwa-samtools-bcftools'

  input:
    file reference from referenceGenomeGz4
    file referenceIndex from referenceIndexes2
    file sample from mergedTrim
  output: file 'reads.sam' into readsUnsorted_trim

  """
  bwa mem -t $THREADS $reference $sample | samtools view -Sbh - -o reads.sam
  """
}

// === SAMTOOLS SORT ===

process sortAlignment {
  container 'biodckr/samtools'

  input: file sam from readsUnsorted
  output: file 'reads.bam' into readsBAM

  """
  samtools sort -@ $THREADS $sam -o reads.bam > reads.bam
  """
}

( readsBAM1,
  readsBAM2 ) = readsBAM.into(2)

process sortAlignment_trim {
  container 'biodckr/samtools'

  input: file sam from readsUnsorted_trim
  output: file 'reads.bam' into readsBAM_trim

  """
  samtools sort $sam -o reads.bam > reads.bam
  """
}

( readsBAM_trim1,
  readsBAM_trim2 ) = readsBAM_trim.into(2)

// === SAMTOOLS INDEX ===

process indexAlignment {
  container 'biodckr/samtools'

  input: file bam from readsBAM1
  output: file 'reads.bam.bai' into readsBAI

  """
  samtools index $bam
  """
}

process indexAlignment_trim {
  container 'biodckr/samtools'

  input: file bam from readsBAM_trim1
  output: file 'reads.bam.bai' into readsBAI_trim

  """
  samtools index $bam
  """
}

// === SAMTOOLS MPILEUP | BCFTOOLS CALL ===

process callVariants {
  container 'bwa-samtools-bcftools'

  input:
    file bam from readsBAM2
    file bai from readsBAI
    file reference from referenceGenomes1
  output: 'variants.vcf'

  """
  samtools mpileup -uf $reference $bam | bcftools call -c - > variants.vcf
  """
}

process callVariants {
  container 'bwa-samtools-bcftools'

  input:
    file bam from readsBAM_trim2
    file bai from readsBAI_trim
    file reference from referenceGenomes2
  output: 'variants.vcf'

  """
  samtools mpileup -uf $reference $bam | bcftools call -c - > variants.vcf
  """
}
