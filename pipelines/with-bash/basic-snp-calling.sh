#!/bin/bash
set -e

# Very basic SNP calling pipeline example for testing/development
# Copyright 2016 Bruno Vieira <mail@bmpvieira.com> - MIT license
# Edited by Julian Mazzitelli <mazzitelli.julian@gmail.com>

# Notes: For simplicity this does not do any quality filtering and uses small
# datasets that will probably not provide any biologically relevant result
# Edit: Now filters using khmer and kmc.

# Requirements: bionode-ncbi, sra-toolkit, bwa, seqtk, samtools, and bcftools
# (tip: check GitHub, Docker, Homebrew or Linuxbrew)

# Workflow parameters
# Salmonella enterica
REFERENCE_NAME='GCA_000006945.2_ASM694v2'
REFERENCE="${REFERENCE_NAME}_genomic.fna.gz"
REFERENCE_URL="http://ftp.ncbi.nlm.nih.gov/genomes/all/$REFERENCE_NAME/$REFERENCE"
READS='2492428'
export READSFQ='ERR1229296'
FILTER_MODE='khmer' # 'khmer'

# Config
export THREADS=2
export MEMORYGB=4
# Need to mkdir this first..
export TMPDIR=$HOME/kmc-temp
TMP=./tmp
export BIN=../bin
ADAPTERS=../adapters


export KMERSIZE=20
export PLOTXMAX=60
export PLOTYMAX=1200000
# Set this based off k-mer plot.
export MINCOVERAGE=8

# === DOWNLOAD ===

# Download reference genome
# This works for some species, others it will download rna_from_genomic
# see: https://github.com/bionode/bionode-ncbi/issues/19
# bionode-ncbi download assembly $SPECIES
# TODO fix rna_from_genomic bug, use bionode-ncbi instead of curl
echo "START reference download"
curl -sO $REFERENCE_URL
echo "END reference download"

# Download sequencing reads
echo "START reads download"
bionode-ncbi download sra $READS > $TMP
echo "END reads download"

# === EXTRACT READS ===
echo "START fastq-dump"
fastq-dump --split-files --skip-technical --gzip $READS/**.sra
echo "END fastq-dump"

# === TRIM READS ADAPTERS ===
# following command is for Illumina paired end reads, so it needs to be adjusted
# for single end or other technologies (e.g., change the adapters/TruSeq3-PE.fa file)
# Adapter files are included with Trimmomatic installation.
# Parameters from quick start guide: http://www.usadellab.org/cms/?page=trimmomatic

echo "START trim reads adapters"
java -jar $BIN/trimmomatic-0.36.jar PE -phred33 \
  $READSFQ\_1.fastq.gz $READSFQ\_2.fastq.gz \
  $READSFQ\_1.trim.pe.fastq.gz $READSFQ\_1.trim.se.fastq.gz \
  $READSFQ\_2.trim.pe.fastq.gz $READSFQ\_2.trim.se.fastq.gz \
  ILLUMINACLIP:$ADAPTERS/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
echo "END trim reads adapters"

# merge paired end fasta files into interleaved fastq
echo "START merge trimmed reads"
seqtk mergepe $READSFQ\_1.trim.pe.fastq.gz $READSFQ\_2.trim.pe.fastq.gz | gzip > $READSFQ.trim.pe.fastq.gz
echo "END merge trimmed reads"

# === FILTER RARE k-mers ===

if [ $FILTER_MODE == 'kmc' ]; then
  if [ ! -d $TMPDIR ]; then
    echo "Need to mkdir $TMPDIR"
    exit 1
  fi
  echo "START filtering with kmc"
  export READS_TO_ALIGN="$READSFQ.trim.pe.kmc.fastq.gz"
  ./filter_kmc.sh
  echo "END filtering with kmc"
elif [ $FILTER_MODE == 'khmer' ]; then
  echo "Using khmer mode"
  export READS_TO_ALIGN="$READSFQ.trim.pe.khmer.fastq.gz"
  ./filter_khmer.sh
  FINAL_OUTPUT="${READS}-khmer.vcf"
else
  echo "No filter mode set. Continuing with reads with adapters trimmed."
  READS_TO_ALIGN="$READSFQ.trim.pe.fastq.gz"
fi

# === READS ALIGNMENT ===

echo "START index reference"
bwa index $REFERENCE
echo "END index reference"

echo "START Align reads to reference genome"
bwa mem -t $THREADS $REFERENCE $READS_TO_ALIGN > $READS.sam
echo "END reads aligned"

# Convert alignment file to binary
echo "START sam -> bam"
samtools view -@ $THREADS -bh $READS.sam  > $READS.unsorted.bam
echo "END sam -> bam"

echo "START Sort alignment file"
samtools sort -@ $THREADS $READS.unsorted.bam -o $READS.bam
echo "START Sort alignment file"

# Note: previous steps could probably also work like:
# samtools view -bh -@ $THREADS $READS.sam | samtools sort -@ $THREADS - $READS.bam
# Or even pipe from bwa step, check: http://www.htslib.org/workflow/

# makes $READS.bam.bai (binary alignment index)
echo "START index alignment"
samtools index $READS.bam
echo "END index alignment"

# === VARIANT CALLING ===

# convert .fna.gz to .fna with bgzip because:
# [fai_load] build FASTA index.
# [fai_build] fail to open the FASTA file 503988/GCA_000315625.1_Guith1_genomic.fna.gz
echo "START decompress reference"
bgzip -@ $THREADS -d $REFERENCE
echo "END decompress reference"

# then multipileup with .fna, and call with bcftools
echo "START call variants"
samtools mpileup -uf `echo $REFERENCE | sed s/\.gz$//` $READS.bam | bcftools call -c - > $FINAL_OUTPUT
echo "END call variants"
