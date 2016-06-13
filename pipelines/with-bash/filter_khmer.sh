#!/bin/bash

# using khmer
# count k-mers
# 25mins
echo "START load-into-counting"
load-into-counting.py -N 4 -k $KMERSIZE -M 8e9 -T $THREADS $READSFQ.trim.pe.fastq.gz.kh  $READSFQ.trim.pe.fastq.gz
echo "END load-into-counting"

# create histogram of counts
echo "START make histogram"
abundance-dist.py $READSFQ.trim.pe.fastq.gz.kh $READSFQ.trim.pe.fastq.gz $READSFQ.trim.pe.fastq.gz.kh.hist
# plot histogram (xmax and ymax might need adjustments so that plot can be interpreted)
# plot-abundance-dist.py --xmax $PLOTXMAX --ymax $PLOTYMAX $READSFQ.trim.pe.fastq.gz.kh.hist $READSFQ.trim.pe.fastq.gz.kh.hist.png
./plotKhmer.R $READSFQ.trim.pe.fastq.gz.kh.hist $READSFQ.trim.pe.fastq.gz.kh.hist.png $PLOTXMAX $PLOTYMAX
echo "END stop histogram"

# INTERACTION REQUIRED: from plot results, pick a minimum frequency cutoff for filter
# EXPECTS THIS TO BE DEFINED FROM PARENT SCRIPT
# MINCOVERAGE=5

# filter low-abundance k-mers from reads (reads are not entirely discarded)
# 14mins
echo "START filter-abound"
filter-abund.py -T $THREADS -C $MINCOVERAGE $READSFQ.trim.pe.fastq.gz.kh -o $READS_TO_ALIGN $READSFQ.trim.pe.fastq.gz
echo "END filter-abound"
