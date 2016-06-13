#!/usr/bin/env Rscript
# Quickly plot KMC histograms with R and ggplot2
# Copyright 2016 Bruno Vieira <mail@bmpvieira.com> - MIT license

# Usage:
#   ./plotKMC.R <input_file> <output_file.png> <x_limit> <y_limit>

library("ggplot2")
args = commandArgs(trailingOnly=TRUE)
print(args[1])
data = read.table(args[1])
names(data) = c("Coverage", "Count")
png(args[2], width=960, height=960)
ggplot(data, aes(x=Coverage, y=Count)) +
  geom_line(colour="blue") +
  coord_cartesian(xlim=c(0, as.integer(args[3])), ylim=c(0, as.numeric(args[4]))) +
  scale_x_continuous(minor_breaks=seq(1,1000,1)) +
  scale_y_continuous(minor_breaks=seq(0,1,0.10))
dev.off()
