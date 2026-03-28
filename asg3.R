# Stats Asg 3

# --------------- Import Libraries ---------------
library(VariantAnnotation)
library(tidyverse)
library(snpStats)
library(cluster)
library(patchwork)

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

gt_all <- bind_rows(fto_eda, tcf_eda, slc_eda)

# AF distribution per population and gene
all_eda_long <- gt_all %>%
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


# --------------- PCA: Preparation ---------------

# Function to prepare genotype matrix
prep_matrix <- function(gt_matrix) {
  # Fill NAs with row mean
  gt_filled <- t(apply(gt_matrix, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  }))
  # Remove monomorphic SNPs
  gt_filtered <- gt_filled[apply(gt_filled, 1, var) > 0, ]
  cat("SNPs remaining after monomorphic filter:", nrow(gt_filtered), "\n")
  return(gt_filtered)
}

fto_mat <- prep_matrix(fto_gt)
tcf_mat <- prep_matrix(tcf_gt)
slc_mat <- prep_matrix(slc_gt)

# Function to run PCA
run_pca <- function(gt_matrix) {
  prcomp(t(gt_matrix), scale. = TRUE)
}

fto_pca <- run_pca(fto_mat)
tcf_pca <- run_pca(tcf_mat)
slc_pca <- run_pca(slc_mat)

# Function to get variance explained
get_var_explained <- function(pca_result) {
  eigenvals <- pca_result$sdev^2
  eigenvals / sum(eigenvals) * 100
}

fto_var <- get_var_explained(fto_pca)
tcf_var <- get_var_explained(tcf_pca)
slc_var <- get_var_explained(slc_pca)

# Scree plot
scree_df <- bind_rows(
  data.frame(PC = 1:20, Variance = fto_var[1:20], gene = "FTO"),
  data.frame(PC = 1:20, Variance = tcf_var[1:20], gene = "TCF7L2"),
  data.frame(PC = 1:20, Variance = slc_var[1:20], gene = "SLC30A8")
)

ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~gene) +
  scale_x_continuous(breaks = 1:20) +
  labs(title = "Scree Plots by Gene",
       x     = "PC Number",
       y     = "Variance Explained (%)") +
  theme_minimal()

# Function to build PCA score dataframes with population labels
build_pca_scores <- function(pca_result, panel, var_explained, gene_name) {
  scores         <- as.data.frame(pca_result$x)
  scores$sample  <- rownames(scores)
  scores         <- left_join(scores, panel, by = "sample")
  scores$gene    <- gene_name
  scores$pc1_var <- round(var_explained[1], 1)
  scores$pc2_var <- round(var_explained[2], 1)
  return(scores)
}

fto_scores <- build_pca_scores(fto_pca, panel, fto_var, "FTO")
tcf_scores <- build_pca_scores(tcf_pca, panel, tcf_var, "TCF7L2")
slc_scores <- build_pca_scores(slc_pca, panel, slc_var, "SLC30A8")

# Quick check
head(fto_scores[, c("sample", "pop", "super_pop")])


# --------------- PCA Plots ---------------
# Superpopulation colours
superpop_colours <- c("EUR" = "coral", "SAS" = "steelblue")

# Subpopulation colours (5 EUR + 5 SAS = 10 subpops)
subpop_colours <- c(
  # EUR subpopulations
  "CEU" = "#e63946", "TSI" = "#f4a261", "FIN" = "#e9c46a",
  "GBR" = "#2a9d8f", "IBS" = "#457b9d",
  # SAS subpopulations  
  "BEB" = "#7b2d8b", "GIH" = "#c77dff", "ITU" = "#9d4edd",
  "PJL" = "#5a189a", "STU" = "#e0aaff"
)

# Plot function for superpopulation
plot_pca_superpop <- function(scores, var_exp, gene_name) {
  ggplot(scores, aes(x = PC1, y = PC2, colour = super_pop)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = superpop_colours) +
    labs(title  = gene_name,
         x      = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y      = paste0("PC2 (", round(var_exp[2], 1), "%)"),
         colour = "Superpopulation") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Plot function for subpopulation (not doing anymore)
plot_pca_subpop <- function(scores, var_exp, gene_name) {
  ggplot(scores, aes(x = PC1, y = PC2, colour = pop)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = subpop_colours) +
    labs(title  = gene_name,
         x      = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y      = paste0("PC2 (", round(var_exp[2], 1), "%)"),
         colour = "Subpopulation") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Plot superpopulation per gene
(plot_pca_superpop(fto_scores, fto_var, "FTO") +
    plot_pca_superpop(tcf_scores, tcf_var, "TCF7L2") +
    plot_pca_superpop(slc_scores, slc_var, "SLC30A8")) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "PCA - EUR vs SAS Superpopulations") &
  theme(legend.position = "bottom")


# Plot subpopulations per gene
(plot_pca_subpop(fto_scores, fto_var, "FTO") +
    plot_pca_subpop(tcf_scores, tcf_var, "TCF7L2") +
    plot_pca_subpop(slc_scores, slc_var, "SLC30A8")) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "PCA - EUR vs SAS Subpopulations")&
  theme(legend.position = "bottom")

# --------------- K-means: Elbow plot ---------------
# Function to calculate within-cluster sum of squares for k=1 to k=10
elbow_plot <- function(pca_result, title) {
  
  # Try k from 1 to 10
  wss <- sapply(1:10, function(k) {
    kmeans(pca_result$x[, 1:10], 
           centers = k, 
           nstart = 25)$tot.withinss # # run 25 times with diff. starting points, keep best. Then extract total within-cluster sum of squares
  })
  data.frame(k = 1:10, wss = wss, gene = title)
}

# Combine genes
elbow_df <- bind_rows(
  elbow_plot(fto_pca, "FTO"),
  elbow_plot(tcf_pca, "TCF7L2"),
  elbow_plot(slc_pca, "SLC30A8")
)

# Plot
ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_point(size = 2) +
  geom_line() +
  facet_wrap(~gene, scales = "free_y") + # Each plot has own axis
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Elbow Plot by Gene",
       x     = "Number of Clusters (k)",
       y     = "Total Within-Cluster SS") +
  theme_minimal()

# --------------- K-means: Calculation ---------------
set.seed(42)

# Run k-means on first 10 PCs
run_kmeans <- function(pca_result, scores, k) {
  km             <- kmeans(pca_result$x[, 1:10], centers = k, nstart = 25)
  scores$cluster <- as.factor(km$cluster)
  return(scores)
}

# Call function for each gene
fto_k <- run_kmeans(fto_pca, fto_scores, 6)
tcf_k <- run_kmeans(tcf_pca, tcf_scores, 9)
slc_k <- run_kmeans(slc_pca, slc_scores, 4)

# Plot function for k-means
plot_kmeans_compare <- function(scores, var_exp, gene_name, k) {
ggplot(scores, aes(PC1, PC2, colour = cluster)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(title  = paste0("K-means of ", gene_name, " gene (k=",k,")"),
       x      = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y      = paste0("PC2 (", round(var_exp[2], 1), "%)"),
       colour = "Cluster") +
  theme_minimal()
}

# Plots per gene
plot_kmeans_compare(fto_k, fto_var, "FTO", 6)
plot_kmeans_compare(tcf_k, tcf_var, "TCF7L2", 9)
plot_kmeans_compare(slc_k, slc_var, "SLC30A8", 4)

# --------------- K-means: Silhouette score ---------------
# Measures how well each point fits its cluster
fto_sil  <- silhouette(as.integer(fto_k$cluster),  dist(fto_pca$x[, 1:10]))
tcf_sil  <- silhouette(as.integer(tcf_k$cluster),  dist(tcf_pca$x[, 1:10]))
slc_sil  <- silhouette(as.integer(slc_k$cluster),  dist(slc_pca$x[, 1:10]))

# Average silhouette width per gene
# Closer to 1 = well clustered, closer to 0 = overlapping clusters
summary(fto_sil)$avg.width
summary(tcf_sil)$avg.width
summary(slc_sil)$avg.width
