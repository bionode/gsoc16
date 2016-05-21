#!/bin/bash
# Very basic SNP calling pipeline example for testing/development
# Copyright 2016 Bruno Vieira <mail@bmpvieira.com> - MIT license

# Notes: For simplicity this does not do any quality filtering and uses small
# datasets that will probably not provide any biologically relevant result

# Requirements: bionode-ncbi, sra-toolkit, bwa, seqtk, samtools, and bcftools
# (tip: check Docker, Homebrew or Linuxbrew)

# Config
THREADS=2
TMPDIR=/tmp
SPECIES=Guillardia theta
REFERENCE=503988/GCA_000315625.1_Guith1_genomic.fna.gz
READS=35526

# Download reference genome
bionode-ncbi download assembly $SPECIES

# Download sequencing reads
bionode-ncbi download sra $READS

# Extract reads
fastq-dump --split-files --skip-technical --gzip $READS/**.sra

# Index reference genome
bwa index $REFERENCE

# Align reads to reference genome
bwa mem -t $THREADS $REFERENCE *.fastq.gz > $READS.sam

# Convert alignment file to binary
samtools view -@ $THREADS -bh $READS.sam  > $READS.unsorted.bam

# Sort alignment file
samtools sort -@ $THREADS $READS.unsorted.bam -o $READS.bam

# Note: previous steps could probably also work like:
# samtools view -bh -@ $THREADS $READS.sam | samtools sort -@ $THREADS - $READS.bam
# Or even pipe from bwa step, check: http://www.htslib.org/workflow/

# Index alignment file
# makes $READS.bam.bai (binary alignment index)
samtools index $READS.bam

# Do variant calling

# convert .fna.gz to .fna with bgzip because:
# [fai_load] build FASTA index.
# [fai_build] fail to open the FASTA file 503988/GCA_000315625.1_Guith1_genomic.fna.gz
bgzip -@ $THREADS -d $REFERENCE

# then multipileup with .fna
samtools mpileup -uf `echo $REFERENCE | sed s/\.gz$//` $READS.bam | bcftools call -c - > $READS.vcf
