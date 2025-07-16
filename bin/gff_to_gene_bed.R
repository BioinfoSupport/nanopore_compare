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


## Select both genes and CDS
gff_gene <- gff[gff$type %in% c("gene", "CDS")]
## If CDS is inside of gene record it will be marked as duplicate and discarded
gff_gene <- gff_gene[!duplicated(gff_gene)]

## Get the nicest possible gene name
gene_names <- rep(NA, length(gff_gene))
if ("Name" %in% colnames(mcols(gff_gene)))
    gene_names[is.na(gene_names)] <- gff_gene$Name[is.na(gene_names)]
if ("gene" %in% colnames(mcols(gff_gene)))
    gene_names[is.na(gene_names)] <- gff_gene$gene[is.na(gene_names)]
if ("locus_tag" %in% colnames(mcols(gff_gene)))
    gene_names[is.na(gene_names)] <- gff_gene$locus_tag[is.na(gene_names)]
if ("ID" %in% colnames(mcols(gff_gene)))
    gene_names[is.na(gene_names)] <- gff_gene$ID[is.na(gene_names)]
gene_names[is.na(gene_names)] <- ""
names(gff_gene) <- gene_names

# Required for the score (column 5) in bed
gff_gene$score <- 0


export.bed(gff_gene, bed_file)
