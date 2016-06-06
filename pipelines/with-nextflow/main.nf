#!/usr/bin/env nextflow

species = [
  'Salmonella-enterica': [
    'referenceURL': 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000988525.2_ASM98852v2/GCA_000988525.2_ASM98852v2_genomic.fna.gz',
    'readsID': '2492428'
  ],
  'Staphylococcus-aureus': [
    'referenceURL': 'http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000013425.1_ASM1342v1/GCA_000013425.1_ASM1342v1_genomic.fna.gz',
    'readsID': '1274026'
  ]
]

echo true

// === DOWNLOAD ===

process downloadReference {
  input: val referenceURL from species.collect { it.value.referenceURL }
  output: file 'reference.genomic.fna.gz' into referenceGenomesZipped1, referenceGenomesZipped2, referenceGenomesZipped3

  """
  curl $referenceURL -o reference.genomic.fna.gz
  """
}

process downloadSRA {
  input: val readsID from species.collect { it.value.readsID }
  output: file '**/*.sra' into reads

  """
  bionode-ncbi download sra $readsID > tmp
  """
}

// === EXTRACT/DECOMPRESS ===

process extractSRA {
  input: file read from reads
  output: file '*.fastq.gz' into samples

  """
  fastq-dump --split-files --skip-technical --gzip $read
  """
}

process decompressReference {
  input: file referenceGenome from referenceGenomesZipped1
  output: file 'reference.genomic.fna' into referenceGenomes

  """
  bgzip -d $referenceGenome --stdout > reference.genomic.fna
  """
}

// === MAPPING ===

process indexReference {
  input: file reference from referenceGenomesZipped2
  output: file '*.gz.*' into referenceIndexes

  """
  bwa index $reference
  """
}

process mem {
  input:
    file reference from referenceGenomesZipped3
    file referenceIndex from referenceIndexes
    file sample from samples
  output: file 'reads.sam' into readsUnsorted

  """
  bwa mem $reference $sample | samtools view -Sbh -o reads.sam
  """
}

process bam {
  input: file sam from readsUnsorted
  output: file 'reads.bam' into readsBAM1, readsBAM2

  """
  samtools sort $sam -o reads.bam
  """
}

process bai {
  input: file bam from readsBAM1
  output: file 'reads.bam.bai' into readsBAI

  """
  samtools index $bam reads.bam.bai
  """
}

process call {
  input:
    file bam from readsBAM2
    file bai from readsBAI
    file reference from referenceGenomes
  output: 'variants.vcf'

  """
  samtools mpileup -uf $reference $bam | bcftools call -c - > variants.vcf
  """
}
