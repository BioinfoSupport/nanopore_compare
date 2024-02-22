#!/bin/env Rscript

library(GenomicAlignments)
library(VariantAnnotation)

annotate_mm2_vranges <- function(cons_vcf_file, cons_bam, min_bam_mapq = 60) {
    bam <- readGAlignments(cons_bam,
                           param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE,
                                                                   isSupplementaryAlignment = FALSE),
                                                what = c("flag", "qual", "mapq")),
                           use.names = TRUE)
    ## bam <- bam[mcols(bam)$mapq >= min_bam_mapq]

    vr <- readVcfAsVRanges(cons_vcf_file)
    
    ## We are haploid, so remove alleles
    ## TODO: Why it was in diploid form anyway?
    vr$GT <- substring(vr$GT,1,1)
    mcols(vr)$GQ <- as(
                 subseq(mcols(bam[vr$QNAME])$qual, vr$QSTART, width=1),
                 "IntegerList")
    mcols(vr)$CONS3Q <- as(
                 subseq(mcols(bam[vr$QNAME])$qual, vr$QSTART, width=3),
                 "IntegerList")
    vr$QUAL <- as.numeric(vr$GQ)

    vr
}
    
annotate_mm2_vcf <- function(cons_vcf_file, cons_bam, min_bam_mapq = 60) {
    v <- readVcf(cons_vcf_file)
    vr <- annotate_mm2_vranges(cons_vcf_file, cons_bam, min_bam_mapq)

    vnew <- asVCF(vr, info=names(info(v)))
    geno(vnew)$AD <- NULL
    geno(vnew)$DP <- NULL
    geno(vnew)$FT <- NULL
    geno(vnew)$QUAL <- NULL
    info(header(vnew)) <- info(header(v))
    geno(header(vnew)) <- geno(header(vnew))["GT", ] # LEave only this in genotype original description
    geno(header(vnew))["GQ","Number"] <- "1"
    geno(header(vnew))["GQ","Type"] <- "Integer"
    geno(header(vnew))["GQ","Description"] <- "Consensus quality score"
    geno(header(vnew))["CONS3Q","Number"] <- "3"
    geno(header(vnew))["CONS3Q","Type"] <- "Integer"
    geno(header(vnew))["CONS3Q","Description"] <- "Consensus quality 3nt"
    vnew
}


ARGS <- commandArgs(trailingOnly = TRUE)

cons_vcf <- ARGS[1]
cons_bam <- ARGS[2]
out <- ARGS[3]

v <- annotate_mm2_vcf(cons_vcf, cons_bam)
writeVcf(v, out)
