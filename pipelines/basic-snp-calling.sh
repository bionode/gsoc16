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

# Toggle comments to change species..

SPECIES='Solenopsis invicta'
REFERENCE=244018/GCA_000188075.1_Si_gnG_genomic.fna.gz
READS=2209130
READSFQ=DRR042176

SPECIES='Guillardia theta'
REFERENCE=503988/GCA_000315625.1_Guith1_genomic.fna.gz
READS=35526

# Salmonella enterica
# SPECIES='GCA_000006945.2'
REFERENCE='GCA_000006945.2_ASM694v2_genomic.fna.gz'
READS='2492428'
READSFQ='ERR1229296'

# Download reference genome
# bionode-ncbi download assembly $SPECIES
# TODO not curl..
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000006945.2_ASM694v2/$REFERENCE

# Download sequencing reads
bionode-ncbi download sra $READS

# Extract reads
fastq-dump --split-files --skip-technical --gzip $READS/**.sra

# Trim reads adapters
# following command is for Illumina paired end reads, so it needs to be adjusted
# for single end or other technologies (e.g., change the adapters/TruSeq3-PE.fa file)
# Adapter files are included with Trimmomatic installation.
# Parameters from quick start guide: http://www.usadellab.org/cms/?page=trimmomatic

# java -jar trimmomatic-0.36.jar PE -threads $THREADS -phred33 \
java -jar trimmomatic-0.36.jar PE -phred33 \
  $READSFQ\_1.fastq.gz $READSFQ\_2.fastq.gz \
  $READSFQ\_1.trim.pe.fastq.gz $READSFQ\_1.trim.se.fastq.gz \
  $READSFQ\_2.trim.pe.fastq.gz $READSFQ\_2.trim.se.fastq.gz \
  ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# For Single end, hardcoded to Guillardia theta
# java -jar trimmomatic-0.36.jar SE -phred33 \
#   SRR070675_2.fastq.gz trim.SRR070675_2.fastq.gz \
#   ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
#   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Filter reads with low-abundance k-mers using khmer
# settings (will need to be adjusted)
KMERSIZE=20
PLOTXMAX=60
PLOTYMAX=1200000
# merge paired end fasta files into interleaved fasta
seqtk mergepe $READSFQ\_1.trim.pe.fastq.gz $READSFQ\_2.trim.pe.fastq.gz | gzip > $READSFQ.trim.pe.fastq.gz
# count k-mers
# 25mins
load-into-counting.py -N 4 -k $KMERSIZE -M 8e9 -T $THREADS $READSFQ.trim.pe.fastq.gz.kh  $READSFQ.trim.pe.fastq.gz
# create histogram of counts
abundance-dist.py $READSFQ.trim.pe.fastq.gz.kh $READSFQ.trim.pe.fastq.gz $READSFQ.trim.pe.fastq.gz.kh.hist
# plot histogram (xmax and ymax might need adjustments so that plot can be interpreted)
# plot-abundance-dist.py --xmax $PLOTXMAX --ymax $PLOTYMAX $READSFQ.trim.pe.fastq.gz.kh.hist $READSFQ.trim.pe.fastq.gz.kh.hist.png
./plotKhmer.R $READSFQ.trim.pe.fastq.gz.kh.hist $READSFQ.trim.pe.fastq.gz.kh.hist.png $PLOTXMAX $PLOTYMAX

# INTERACTION REQUIRED: from plot results, pick a minimum frequency cutoff for filter
MINCOVERAGE=5
# filter low-abundance k-mers from reads (reads are not entirely discarded)
# 14mins
filter-abund.py -T $THREADS -C $MINCOVERAGE $READSFQ.trim.pe.fastq.gz.kh -o $READSFQ.trim.pe.khmer.fastq.gz $READSFQ.trim.pe.fastq.gz

# Faster alternative k-mer filtering using refresh-bio/KMC
# this methods discards reads so in a paired end dataset it might break some pairs
# these orphan reads might need to be removed otherwise this could cause issues downstream
# count k-mers
kmc -k$KMERSIZE -m$MEMORYGB -t$THREADS $READSFQ.trim.pe.fastq.gz $READSFQ.trim.pe.kmc .
# make histogram
kmc_tools histogram $READSFQ.trim.pe.kmc $READSFQ.trim.pe.kmc.hist.txt
# plot histogram (R script not included with kmc)
plotKMC.R $READSFQ.trim.pe.kmc.hist.txt $READSFQ.trim.pe.kmc.hist.png $PLOTXMAX $PLOTYMAX
# INTERACTION REQUIRED: from plots results, pick a minimum frequency cutoff for filter
MINCOVERAGE=8
# Filter rare k-mers in one step (slower but uses less disk space)
kmc_tools filter $READSFQ.trim.pe.kmc -cx$MINCOVERAGE $READSFQ.trim.pe.fastq.gz -ci0 -cx0 $READSFQ.trim.pe.kmc.fastq.gz
# Or, filter in two steps (faster but requires more disk space)
kmc_tools reduce $READSFQ.trim.pe.kmc  -cx$MINCOVERAGE $READSFQ.rare_k-mers.kmc
kmc_tools filter $READSFQ.rare_k-mers.kmc $READSFQ.trim.pe.fastq.gz -ci0 -cx0 $READSFQ.trim.pe.kmc.fastq.gz
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
