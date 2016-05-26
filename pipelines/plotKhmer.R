#!/usr/bin/env Rscript
# Plot khmer histogram

args <- commandArgs(trailingOnly=TRUE)

fileName <- args[1]
output <- args[2]
maxx <- as.integer(args[3])
maxy <- as.integer(args[4])

data <- read.table(fileName, header=TRUE, sep=",")


png(output, width=960, height=960)
plot(data[,1], data[,2], xlim=c(0,maxx), ylim=c(0, maxy), xlab=colnames(data)[1], ylab=colnames(data)[2])
dev.off()
