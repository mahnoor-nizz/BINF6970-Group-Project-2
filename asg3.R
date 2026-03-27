# Stats Asg 3
library(VariantAnnotation)
library(tidyverse)
library(snpStats)

#BiocManager::install("snpStats")

# Load files
fto_vcf <- readVcf("fto_chr16.vcf.gz")
tcf_vcf <- readVcf("tcf_chr10.vcf.gz")
slc_vcf <- readVcf("slc_chr8.vcf.gz")

# Check number of rows
nrow(fto_vcf)
nrow(tcf_vcf)
nrow(slc_vcf)

# Check number of columns
ncol(geno(fto_vcf)$GT)
ncol(geno(tcf_vcf)$GT)
ncol(geno(slc_vcf)$GT)

# Check the header info fields
info(header(fto_vcf))

# Check if AF columns exist like EAS_AF, EUR_AF, SAS_AF etc.
names(info(fto_vcf))

# Convert into dataframes
fto_info <- as.data.frame(info(fto_vcf))
tcf_info <- as.data.frame(info(tcf_vcf))
slc_info <- as.data.frame(info(slc_vcf))

# Splits the multi alleles based on the higher frequency
fto_keep <- sapply(fto_info$SAS_AF, max) > 0.001 & 
  sapply(fto_info$EUR_AF, max) > 0.001

tcf_keep <- sapply(tcf_info$SAS_AF, max) >= 0.001 & 
  sapply(tcf_info$EUR_AF, max) >= 0.001

slc_keep <- sapply(slc_info$SAS_AF, max) > 0.001 & 
  sapply(slc_info$EUR_AF, max) > 0.001

desired_cols <- c( "AN", "AC", "AF", "SAS_AF", "EUR_AF", "VT", "DP", "EX_TARGET")

# Subset the dataframe
fto_filtered <- fto_info[fto_keep, desired_cols]
tcf_filtered <- tcf_info[tcf_keep, desired_cols]
slc_filtered <- slc_info[slc_keep, desired_cols]

head(tcf_filtered)

# Download the panel file
# download.file(
#   "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
#   destfile = "integrated_call_samples_v3.20130502.ALL.panel"
# )

# Load it
panel <- read.table("integrated_call_samples_v3.20130502.ALL.panel",
                    header = TRUE,
                    col.names = c("sample", "pop", "super_pop", "gender"))

# Quick check
table(panel$super_pop)

# Subset to EUR and SAS
eur_sas_samples <- panel$sample[panel$super_pop %in% c("EUR", "SAS")]
length(eur_sas_samples)

# Sanity check it matches the VCF
head(colnames(fto_vcf))

sum(colnames(fto_vcf) %in% eur_sas_samples)

# Get the genotype for the LD
fto_vcf_sub <- fto_vcf[fto_keep, ]
tcf_vcf_sub <- tcf_vcf[tcf_keep, ]
slc_vcf_sub <- slc_vcf[slc_keep, ]

geno(fto_vcf)

# Function to convert the genotype to numbers
vcf_to_numeric <- function(vcf) {
  gt <- geno(vcf)$GT
  gt <- gsub("\\|", "/", gt)
  
  gt_num <- matrix(NA_real_, nrow=nrow(gt), ncol=ncol(gt),
                   dimnames=dimnames(gt))
  
  gt_num[gt == "0/0"] <- 0
  gt_num[gt %in% c("0/1", "1/0")] <- 1
  gt_num[gt == "1/1"] <- 2
  
  return(gt_num)
}

vcf_to_numeric <- function(vcf) {
  gt <- geno(vcf)$GT
  gt <- gsub("\\|", "/", gt)  # handle phased genotypes
  
  gt_num <- matrix(NA_real_, nrow=nrow(gt), ncol=ncol(gt),
                   dimnames=dimnames(gt))
  
  gt_num[gt == "0/0"] <- 0
  gt_num[gt %in% c("0/1", "1/0")] <- 1
  gt_num[gt == "1/1"] <- 2
  
  return(gt_num)
}

fto_gt <- vcf_to_numeric(fto_vcf_sub)
tcf_gt <- vcf_to_numeric(tcf_vcf_sub)
slc_gt <- vcf_to_numeric(slc_vcf_sub)

fto_gt
tcf_gt
slc_gt

# Convert to SNPmatrix
fto_snp <- as(fto_gt, "SnpMatrix")
tcf_snp <- as(tcf_gt, "SnpMatrix")
slc_snp <- as(slc_gt, "SnpMatrix")

fto_ld <- ld(fto_snp, depth = ncol(fto_snp), stats = "R.squared")
tcf_ld <- ld(tcf_snp, depth = ncol(tcf_snp), stats = "R.squared")
slc_ld <- ld(slc_snp, depth = ncol(slc_snp), stats = "R.squared")

fto_ld

summary(as.vector(fto_ld))
max(fto_ld, na.rm = TRUE)

ld_df <- as.data.frame(as.table(as.matrix(fto_ld)))
colnames(ld_df) <- c("SNP1", "SNP2", "R2")

ggplot(ld_df, aes(SNP1, SNP2, fill = R2)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text = element_blank()) +
  labs(title = "FTO LD")
