#!/bin/bash
# Very basic SNP calling pipeline example for testing/development
# Copyright 2016 Bruno Vieira <mail@bmpvieira.com> - MIT license

# Notes: For simplicity this does not do any quality filtering and uses small
# datasets that will probably not provide any biologically relevant result

# Requirements: bionode-ncbi, sra-toolkit, bwa, seqtk, samtools, and bcftools
# (tip: check Docker, Homebrew or Linuxbrew)

# Config
THREADS=2
MEMORYGB=4
TMPDIR=/tmp
SPECIES='Guillardia theta'
REFERENCE=503988/GCA_000315625.1_Guith1_genomic.fna.gz
READS=35526

# Download reference genome
bionode-ncbi download assembly $SPECIES

# Download sequencing reads
bionode-ncbi download sra $READS

# Extract reads
fastq-dump --split-files --skip-technical --gzip $READS/**.sra

# Trim reads adapters
# following command is for Illumina paired end reads, so it needs to be adjusted
# for single end or other technologies (e.g., change the adapters/TruSeq3-PE.fa file)
# Adapter files are included with Trimmomatic installation.
# Parameters from quick start guide: http://www.usadellab.org/cms/?page=trimmomatic
trimmomatic PE -threads $THREADS -phred33 \
  $READS\_1.fastq.gz $READS\_2.fastq.gz \
  $READS\_1.trim.pe.fastq.gz $READS\_1.trim.se.fastq.gz \
  $READS\_2.trim.pe.fastq.gz $READS\_2.trim.se.fastq.gz \
  ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Filter reads with low-abundance k-mers using khmer
# settings (will need to be adjusted)
KMERSIZE=20
PLOTXMAX=60
PLOTYMAX=12000000
# merge paired end fasta files into interleaved fasta
seqtk mergepe $READS\_1.trim.pe.fastq.gz $READS\_2.trim.pe.fastq.gz | gzip > $READS.trim.pe.fastq.gz
# count k-mers
load-into-counting.py -N 4 -k $KMERSIZE -x 8e9 -T $THREADS $READS.trim.pe.fastq.gz.kh  $READS.trim.pe.fastq.gz
# create histogram of counts
abundance-dist.py $READS.trim.pe.fastq.gz.kh $READS.trim.pe.fastq.gz $READS.trim.pe.fastq.gz.kh.hist
# plot histogram (xmax and ymax might need adjustments so that plot can be interpreted)
plot-abundance-dist.py --xmax $PLOTXMAX --ymax $PLOTYMAX $READS.trim.pe.fastq.gz.kh.hist $READS.trim.pe.fastq.gz.kh.hist.png
# INTERACTION REQUIRED: from plot results, pick a minimum frequency cutoff for filter
MINCOVERAGE=8
# filter low-abundance k-mers from reads (reads are not entirely discarded)
filter-abund.py -T $THREADS -C $MINCOVERAGE $READS.trim.pe.fastq.gz.kh -o $READS.trim.pe.khmer.fastq.gz $READS.trim.pe.fastq.gz

# Faster alternative k-mer filtering using refresh-bio/KMC
# this methods discards reads so in a paired end dataset it might break some pairs
# these orphan reads might need to be removed otherwise this could cause issues downstream
# count k-mers
kmc -k$KMERSIZE -m$MEMORYGB -t$THREADS $READS.trim.pe.fastq.gz $READS.trim.pe.kmc .
# make histogram
kmc_tools histogram $READS.trim.pe.kmc $READS.trim.pe.kmc.hist.txt
# plot histogram (R script not included with kmc)
plotKMC.R $READS.trim.pe.kmc.hist.txt $READS.trim.pe.kmc.hist.png $PLOTXMAX $PLOTYMAX
# INTERACTION REQUIRED: from plots results, pick a minimum frequency cutoff for filter
MINCOVERAGE=8
# Filter rare k-mers in one step (slower but uses less disk space)
kmc_tools filter $READS.trim.pe.kmc -cx$MINCOVERAGE $READS.trim.pe.fastq.gz -ci0 -cx0 $READS.trim.pe.kmc.fastq.gz
# Or, filter in two steps (faster but requires more disk space)
kmc_tools reduce $READS.trim.pe.kmc  -cx$MINCOVERAGE $READS.rare_k-mers.kmc
kmc_tools filter $READS.rare_k-mers.kmc $READS.trim.pe.fastq.gz -ci0 -cx0 $READS.trim.pe.kmc.fastq.gz
# The previous steps might require some explanation. We are filtering from the kmc k-mers table all k-mers
# above our limit with the first -cx (so keeping the rare ones, this is counter intuitive). Then, we remove
# from the fastq file all the reads that have a k-mer present in the filtered table, i.e., a rare k-mer (the -ci0 and -cx0).

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
