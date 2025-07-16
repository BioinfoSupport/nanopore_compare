#!/bin/env Rscript
###
### A simple gff to bed converter, extracting full genes
###
###

ARGS <- commandArgs(trailingOnly = TRUE)
if (length(ARGS) < 2) {
    stop("Usage: gff_to_gene_bed.R source.gff target.bed")
}

gff_file <- ARGS[1]
bed_file <- ARGS[2]

suppressPackageStartupMessages({
    library(rtracklayer)
})


gff <- import.gff(gff_file)


gff_gene <- gff[gff$type == "gene"]

names(gff_gene) <- ifelse(is.na(gff_gene$gene), gff_gene$Name, gff_gene$gene)

# Required for the score (column 5) in bed
gff_gene$score <- 0


export.bed(gff_gene, bed_file)
