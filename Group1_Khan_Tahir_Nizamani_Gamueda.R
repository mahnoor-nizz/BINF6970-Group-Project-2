### =================
# BINF6970 Assignment 3
# Analysis of genotype data from 1000 Genomes Project

### === PACKAGES USED ========
library(VariantAnnotation)
library(TVTB)
library(HardyWeinberg)
library(dplyr)
library(ggplot2)


### === DATA INPUT AND PREPROCESSING ========
# Input VCF files
vcf_fto <- readVcf("fto_chr16.vcf.gz")
vcf_tcf <- readVcf("tcf_chr10.vcf.gz")
vcf_slc <- readVcf("slc_chr8.vcf.gz")

# Input panel files
panel <- read.table("integrated_call_samples_v3.20130502.ALL.panel", header = T)
panel <- panel[panel$super_pop %in% c("SAS", "EUR"), ]
table(panel$pop)

# Check number of SNPs pre-filters
info(vcf_fto) # variants = 11355
info(vcf_tcf) # variants = 7641
info(vcf_slc) # variants = 5873

# Check multiallelic sites
table(elementNROWS(alt(vcf_fto))) # 56 multiallelic sites
table(elementNROWS(alt(vcf_tcf))) # 51 multiallelic sites
table(elementNROWS(alt(vcf_slc))) # 16 multiallelic sites

# Create filter rules (filter out rare alleles)
info_rules <- VcfInfoRules(list(
  af_filt = function(vcf) {
    sas_af <- sapply(vcf$SAS_AF, max)
    eur_af <- sapply(vcf$EUR_AF, max)
    sas_af > 0.001 & eur_af > 0.001
  }
))

# Filtering out rare alleles and checking number of variants post-filter
vcf_fto_filt <- subsetByFilter(vcf_fto, info_rules)
vcf_fto_filt # 1696 variants

vcf_tcf_filt <- subsetByFilter(vcf_tcf, info_rules)
vcf_tcf_filt # 1161 variants

vcf_slc_filt <- subsetByFilter(vcf_slc, info_rules)
vcf_slc_filt # 534 variants

# Function to annotate subpopulations
annotate <- function(vcf) {
  gt_matrix <- geno(vcf)$GT
  samples <- intersect(colnames(gt_matrix), panel$sample)
  
  gt_matrix_sub <- gt_matrix[, samples]
  panel_match <- panel[match(colnames(gt_matrix_sub), panel$sample), ]
  
  sample_meta <- data.frame(sample = panel_match$sample,
                            pop = panel_match$pop,
                            super_pop = panel_match$super_pop,
                            row.names = panel_match$sample)
  return(list(gt_matrix = gt_matrix_sub,
              sample_meta = sample_meta))
}

# Annotated genotype matrices
result_fto <- annotate(vcf_fto_filt)
gt_fto <- result_fto$gt_matrix

result_tcf <- annotate(vcf_tcf_filt)
gt_tcf <- result_tcf$gt_matrix

result_slc <- annotate(vcf_slc_filt)
gt_slc <- result_slc$gt_matrix

# Combining genotype matrices into one 
gt_all <- rbind(gt_fto, gt_tcf, gt_slc)
dim(gt_all) # 3391 SNPs across 3 genes

# Convert genotype to numeric (made function because I don't know if we'll be using separate matrices later)
numeric_convert <- function(gt) {
  gt <- gsub("\\|", "/", gt)
  
  gt_num <- matrix(NA_real_, 
               nrow = nrow(gt),
               ncol = ncol(gt),
               dimnames = dimnames(gt))
  
  gt_num[gt == "0/0"] <- 0
  gt_num[gt %in% c("0/1", "1/0")] <- 1
  gt_num[gt == "1/1"] <- 2
  
  return(gt_num)
}

gt_all_num <- numeric_convert(gt_all)

# Check if conversion worked
gt_all[1:5, 1:5]
gt_all_num[1:5, 1:5]



### === Exploratory Data Analysis ========
# Exploratory PCA
# Account for N/a values in matrix
gt_all_num <- t(apply(gt_all_num, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  x
}))
sum(is.na(gt_all_num)) # all N/A's converted

# Filter monomoprhic SNPs
gt_all_num <- gt_all_num[apply(gt_all_num, 1,var) > 0, ]
nrow(gt_all_num) # 3391 SNPs to 3307 SNPs

# Run PCA
pca_gt_all <- prcomp(t(gt_all_num),
                     scale. = T)
summary(pca_gt_all)

# Calculate variance explained
eigenvals <- pca_gt_all$sdev^2
var_explained <- eigenvals / sum(eigenvals) * 100

# For YUMNA: idk if its variance explained or eigenvalues
scree_df <- data.frame(
  PC = 1:20,
  Variance = var_explained[1:20]
)

ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_continuous(breaks = 1:20) +
  labs(title = "Scree Plot",
       x = "PC Number",
       y = "Variance Explained (%)") +
  theme_minimal()

# Extract PCA scores
pca_scores <- as.data.frame(pca_gt_all$x)
pca_scores$sample <- panel$sample
pca_scores$pop <- panel$pop
pca_scores$super_pop <- panel$super_pop


# PCA plot for superpopulation
ggplot(pca_scores, aes(x = PC1, y = PC2, colour = super_pop)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "PCA of Superpopulation (SAS vs EUR)",
       x = paste0("PC1 (", round(var_explained[1], 1), "% variance explained)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "% variance explained)"),
       colour = "Superpopulation") +
  theme_minimal()
# some overlap, but there are some distinct clusters

# PCA plots for subpopulation
pca_EUR <- pca_scores[pca_scores$super_pop == "EUR", ]
pca_SAS <- pca_scores[pca_scores$super_pop == "SAS", ]

ggplot(pca_SAS, aes(x = PC1, y = PC2, colour = pop)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title  = "PCA of EUR Subpopulations",
       x = paste0("PC1 (", round(var_explained[1], 1), "% variance explained)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "% variance explained)"),
       colour = "Subpopulation") +
  theme_minimal()

ggplot(pca_EUR, aes(x = PC1, y = PC2, colour = pop)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title  = "PCA of EUR Subpopulations",
       x = paste0("PC1 (", round(var_explained[1], 1), "% variance explained)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "% variance explained)"),
       colour = "Subpopulation") +
  theme_minimal()
# large overlap and no distinct clustering

# Hardy Weinberg Equilibrium Test



### === CLUSTERING ANALYSIS ========


### === CLASSIFICATION ANALYSIS ========
