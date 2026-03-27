### =================
# BINF6970 Assignment 3
# Analysis of genotype data from 1000 Genomes Project

### === PACKAGES USED ========
library(VariantAnnotation)
library(TVTB)
library(HardyWeinberg)
library(dplyr)


### === DATA INPUT AND PREPROCESSING ========
vcf_fto <- readVcf("fto_chr16.vcf.gz")
vcf_tcf <- readVcf("tcf_chr10.vcf.gz")
vcf_slc <- readVcf("slc_chr8.vcf.gz")

info(vcf_fto) # variants = 11355
info(vcf_tcf) # variants = 7641
info(vcf_slc) # variants = 5873

table(elementNROWS(alt(vcf_fto))) # 56 multiallelic sites
table(elementNROWS(alt(vcf_tcf))) # 51 multiallelic sites
table(elementNROWS(alt(vcf_slc))) # 16 multiallelic sites

### === EDA ========
info_rules <- VcfInfoRules(list(
  af = function(vcf) {
    af <- sapply(vcf$AF, max)
    af > 0.001 & af < 0.99 # adjust filter for monomorphic alleles if needed
  }
))

vcf_fto_filt <- subsetByFilter(vcf_fto, info_rules)
vcf_fto_filt # 4116 variants

vcf_tcf_filt <- subsetByFilter(vcf_tcf, info_rules)
vcf_tcf_filt # 2796 variants

vcf_slc_filt <- subsetByFilter(vcf_slc, info_rules)
vcf_slc_filt # 1658 variants


# Hardy Weinberg Equilibrium Test



### === CLUSTERING ANALYSIS ========


### === CLASSIFICATION ANALYSIS ========
