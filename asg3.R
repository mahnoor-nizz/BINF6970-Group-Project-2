# Stats Asg 3

# --------------- Import Libraries ---------------
library(VariantAnnotation)
library(tidyverse)
library(snpStats)

#BiocManager::install("snpStats")


# --------------- Load Files and Explore ---------------
# Load files
vcf_fto <- readVcf("fto_chr16.vcf.gz")
vcf_tcf <- readVcf("tcf_chr10.vcf.gz")
vcf_slc <- readVcf("slc_chr8.vcf.gz")

# Check number of rows
nrow(vcf_fto)
nrow(vcf_tcf)
nrow(vcf_slc)

# Check number of columns
ncol(geno(vcf_fto)$GT)
ncol(geno(vcf_tcf)$GT)
ncol(geno(vcf_slc)$GT)

# Check the header info fields
info(header(vcf_fto))

# Check if AF columns exist like EAS_AF, EUR_AF, SAS_AF etc.
names(info(vcf_fto))

# Load panel file
panel <- read.table("integrated_call_samples_v3.20130502.ALL.panel",
                    header = TRUE,
                    col.names = c("sample", "pop", "super_pop", "gender"))

# Quick check
table(panel$super_pop)

# --------------- Filter Data ---------------

# Convert into dataframes
fto_info <- as.data.frame(info(vcf_fto))
tcf_info <- as.data.frame(info(vcf_tcf))
slc_info <- as.data.frame(info(vcf_slc))

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

# Subset panel file to EUR and SAS
eur_sas_samples <- panel$sample[panel$super_pop %in% c("EUR", "SAS")]
length(eur_sas_samples)

# Subset variants and samples
fto_vcf_sub <- vcf_fto[fto_keep, colnames(vcf_fto) %in% eur_sas_samples]
tcf_vcf_sub <- vcf_tcf[tcf_keep, colnames(vcf_tcf) %in% eur_sas_samples]
slc_vcf_sub <- vcf_slc[slc_keep, colnames(vcf_slc) %in% eur_sas_samples]

# --------------- EDA: Allele Frequency Histograms ---------------
# Extract SAS and EUR AF for each gene into dataframes
make_eda_df <- function(info_df, keep_idx, gene_name) {
  data.frame(
    SAS_AF = sapply(info_df[keep_idx, "SAS_AF"], "[", 1),
    EUR_AF = sapply(info_df[keep_idx, "EUR_AF"], "[", 1),
    gene   = gene_name
  )
}

fto_eda <- make_eda_df(fto_info, fto_keep, "FTO")
tcf_eda <- make_eda_df(tcf_info, tcf_keep, "TCF7L2")
slc_eda <- make_eda_df(slc_info, slc_keep, "SLC30A8")

all_eda <- bind_rows(fto_eda, tcf_eda, slc_eda)

# AF distribution per population and gene
all_eda_long <- all_eda %>%
  pivot_longer(cols = c(SAS_AF, EUR_AF),
               names_to  = "population",
               values_to = "AF") %>%
  mutate(population = recode(population,
                             "SAS_AF" = "SAS",
                             "EUR_AF" = "EUR"))

# Plot the AF of each population per gene
ggplot(all_eda_long, aes(x = AF, fill = population)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~gene) +
  scale_fill_manual(values = c("SAS" = "steelblue", "EUR" = "coral")) +
  labs(title = "Allele Frequency Distribution",
       x = "Allele Frequency", y = "Count")


# --------------- EDA: Calculate LD ---------------

# Function to convert the genotype to numbers
vcf_to_numeric <- function(vcf) {
  # Get genotype
  gt <- geno(vcf)$GT
  
  # Replace | with /
  gt <- gsub("\\|", "/", gt)
  
  # Get empty matrix
  gt_num <- matrix(NA_real_, nrow=nrow(gt), ncol=ncol(gt),
                   dimnames=dimnames(gt))
  
  # Get alleles
  gt_num[gt == "0/0"] <- 0
  gt_num[gt %in% c("0/1", "1/0")] <- 1
  gt_num[gt == "1/1"] <- 2
  
  return(gt_num)
}

# Call function for each gene
fto_gt <- vcf_to_numeric(fto_vcf_sub)
tcf_gt <- vcf_to_numeric(tcf_vcf_sub)
slc_gt <- vcf_to_numeric(slc_vcf_sub)

# Convert to SNPmatrix
fto_snp <- as(fto_gt, "SnpMatrix")
tcf_snp <- as(tcf_gt, "SnpMatrix")
slc_snp <- as(slc_gt, "SnpMatrix")

# Compute R squared
fto_ld <- ld(fto_snp, depth = ncol(fto_snp), stats = "R.squared")
tcf_ld <- ld(tcf_snp, depth = ncol(tcf_snp), stats = "R.squared")
slc_ld <- ld(slc_snp, depth = ncol(slc_snp), stats = "R.squared")

# Summary of R squared value
summary(as.vector(fto_ld))
max(fto_ld, na.rm = TRUE)

# Create a df for the heatmap
ld_df <- as.data.frame(as.table(as.matrix(fto_ld)))
colnames(ld_df) <- c("SNP1", "SNP2", "R2")

# Plot
ggplot(ld_df, aes(SNP1, SNP2, fill = R2)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text = element_blank()) +
  labs(title = "FTO LD")

# --------------- K means ---------------