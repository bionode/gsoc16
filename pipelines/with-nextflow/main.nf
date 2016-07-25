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
TMPDIR = 'kmc-temp'

KMERSIZE = 20
PLOTXMAX = 60
PLOTYMAX = 1200000
MINCOVERAGE = 5

THREADS = 20
MEMORYGB = 4

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
referenceGenomeGz.into { referenceGenomeGz1;
                         referenceGenomeGz2;
                         referenceGenomeGz3;
                         referenceGenomeGz4 }

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

samples.into { samples1; samples2 }

process decompressReference {
  container 'biodckrdev/htslib'

  input: file referenceGenome from referenceGenomeGz1
  output: file 'reference.genomic.fna' into referenceGenomes

  """
  bgzip -d $referenceGenome --stdout > reference.genomic.fna
  """
}

referenceGenomes.into { referenceGenomes1; referenceGenomes2 }

// === TRIMMING ===

process trim {
  // container 'thejmazz/polyglot-ngs-01'

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

mergedTrim.into { mergedTrim1; mergedTrim2 }

// === FILTERING ===
process filterKMC {
  // container 'thejmazz/polyglot-ngs-01'

  input: file reads from mergedTrim1
  output: file 'reads.trim.pe.kmc.fastq.gz' into reads_kmc

  // skipping histogram for now

  """
  mkdir $TMPDIR
  kmc -k${KMERSIZE} -m${MEMORYGB} -t${THREADS} $reads reads.trim.pe.kmc ${TMPDIR}
  kmc_tools filter reads.trim.pe.kmc -cx${MINCOVERAGE} $reads -ci0 -cx0 reads.trim.pe.kmc.fastq.gz
  """
}

process filterKHMER {
  // container 'thejmazz/polyglot-ngs-01'

  input: file reads from mergedTrim2
  output: file 'reads.trim.pe.khmer.fastq.gz' into reads_khmer

  // skipping histogram plot

  """
  load-into-counting.py -N 4 -k ${KMERSIZE} -M ${MEMORYGB}e9 -T ${THREADS} reads.trim.pe.fastq.gz.kh $reads
  abundance-dist.py reads.trim.pe.fastq.gz.kh $reads reads.trim.pe.fastq.gz.kh.hist
  filter-abund.py -T ${THREADS} -C ${MINCOVERAGE} reads.trim.pe.fastq.gz.kh -o reads.trim.pe.khmer.fastq.gz $reads
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

referenceIndexes.into { referenceIndexes1; referenceIndexes2 }

// === BWA MEM ===

process alignReads_kmc {
  container 'bwa-samtools-bcftools'

  input:
    file reference from referenceGenomeGz3
    file referenceIndex from referenceIndexes1
    file sample from reads_kmc
  output: file 'reads.sam' into readsUnsorted_kmc

  """
  bwa mem -t $THREADS $reference $sample | samtools view -Sbh - -o reads.sam
  """
}

process alignReads_khmer {
  container 'bwa-samtools-bcftools'

  input:
    file reference from referenceGenomeGz4
    file referenceIndex from referenceIndexes2
    file sample from reads_khmer
  output: file 'reads.sam' into readsUnsorted_khmer

  """
  bwa mem -t $THREADS $reference $sample | samtools view -Sbh - -o reads.sam
  """
}

// === SAMTOOLS SORT ===

process sortAlignment_kmc {
  container 'biodckr/samtools'

  input: file sam from readsUnsorted_kmc
  output: file 'reads.bam' into readsBAM_kmc

  """
  samtools sort $sam -o reads.bam > reads.bam
  """
}

readsBAM_kmc.into { readsBAM_kmc1; readsBAM_kmc2 }

process sortAlignment_khmer {
  container 'biodckr/samtools'

  input: file sam from readsUnsorted_khmer
  output: file 'reads.bam' into readsBAM_khmer

  """
  samtools sort $sam -o reads.bam > reads.bam
  """
}

readsBAM_khmer.into {readsBAM_khmer1; readsBAM_khmer2 }

// === SAMTOOLS INDEX ===

process indexAlignment_kmc {
  container 'biodckr/samtools'

  input: file bam from readsBAM_kmc1
  output: file 'reads.bam.bai' into readsBAI_kmc

  """
  samtools index $bam
  """
}

process indexAlignment_khmer {
  container 'biodckr/samtools'

  input: file bam from readsBAM_khmer1
  output: file 'reads.bam.bai' into readsBAI_khmer

  """
  samtools index $bam
  """
}

// === SAMTOOLS MPILEUP | BCFTOOLS CALL ===

process callVariants_kmc {
  container 'bwa-samtools-bcftools'

  input:
    file bam from readsBAM_kmc2
    file bai from readsBAI_kmc
    file reference from referenceGenomes1
  output: 'variants.vcf'

  """
  samtools mpileup -uf $reference $bam | bcftools call -c - > variants.vcf
  """
}

process callVariants_khmer {
  container 'bwa-samtools-bcftools'

  input:
    file bam from readsBAM_khmer2
    file bai from readsBAI_khmer
    file reference from referenceGenomes2
  output: 'variants.vcf'

  """
  samtools mpileup -uf $reference $bam | bcftools call -c - > variants.vcf
  """
}
