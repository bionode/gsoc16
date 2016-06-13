#!/bin/bash

# Faster alternative k-mer filtering using refresh-bio/KMC
# this methods discards reads so in a paired end dataset it might break some pairs
# these orphan reads might need to be removed otherwise this could cause issues downstream

# count k-mers
# NOTE need to increase limit for # of opened files on osx: ulimit -n 2048
echo "START count k-mers"
kmc -k$KMERSIZE -m$MEMORYGB -t$THREADS $READSFQ.trim.pe.fastq.gz $READSFQ.trim.pe.kmc $TMPDIR
echo "END count k-mers"

echo "START make histogram"
kmc_tools histogram $READSFQ.trim.pe.kmc $READSFQ.trim.pe.kmc.hist.txt
# plot histogram (R script not included with kmc)
$BIN/plotKMC.R $READSFQ.trim.pe.kmc.hist.txt $READSFQ.trim.pe.kmc.hist.png $PLOTXMAX $PLOTYMAX
echo "END plotted histogram"

# INTERACTION REQUIRED: from plots results, pick a minimum frequency cutoff for filter
# Set this based off of k-mer plot. For now, we will hardcode 8, which was
# observed previously.
# export MINCOVERAGE=8

if [ -z ${MINCOVERAGE} ]; then
  echo "Need to set MINCOVERAGE"
  exit 1
fi

# Filter rare k-mers in one step (slower but uses less disk space)
echo "START filtering rare k-mers"
kmc_tools filter $READSFQ.trim.pe.kmc -cx$MINCOVERAGE $READSFQ.trim.pe.fastq.gz -ci0 -cx0 $READSFQ.trim.pe.kmc.fastq.gz
echo "END filtering rare k-mers"
# Or, filter in two steps (faster but requires more disk space)
# kmc_tools reduce $READSFQ.trim.pe.kmc  -cx$MINCOVERAGE $READSFQ.rare_k-mers.kmc
# kmc_tools filter $READSFQ.rare_k-mers.kmc $READSFQ.trim.pe.fastq.gz -ci0 -cx0 $READSFQ.trim.pe.kmc.fastq.gz

# The previous steps might require some explanation. We are filtering from the kmc k-mers table all k-mers
# above our limit with the first -cx (so keeping the rare ones, this is counter intuitive). Then, we remove
# from the fastq file all the reads that have a k-mer present in the filtered table, i.e., a rare k-mer (the -ci0 and -cx0).
