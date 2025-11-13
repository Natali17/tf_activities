### PART 1. Setup and data loading
library(dplyr)
library(limma) # For batch correction
library(edgeR) # For DGEList, TMM normalization and CPM calculation
library(tibble)
library(ggplot2) # For visualization

# Set the working directory
setwd("//wsl.localhost/Ubuntu/home/natal/maradoner_project")

# Load tumor and cell lines expression metadata
tumor_clusters <- read.csv("C:/Users/natal/PycharmProjects/contests/promoters/october2025/TCGA_only_labeled_meta.csv", header = TRUE, row.names = 1)
cell_lines_clusters <- read.csv("C:/Users/natal/PycharmProjects/contests/promoters/october2025/Cell_lines_protocols.csv", header = TRUE, row.names = 1)

# Load TCGA and cell lines expression data
tcga <- read.table("C:/users/natal/PycharmProjects/contests/promoters/october2025/tcga_total_counts.tsv", header = 1, row.names = 1)
cell_lines <- read.table("C:/users/natal/PycharmProjects/contests/promoters/october2025/cell_lines_total_counts.tsv", header = 1, row.names = 1)

### PART 2. Data filtering and merging
# Define gene types to retain for analysis
good_types <- c('IG_C_gene','IG_V_gene','protein_coding','TR_V_gene','transcribed_unitary_pseudogene')

# --- Filter and summarize TCGA counts ---
# 1. Remove mitochondrial genes (MT-), 2. Filter by good gene type, 3. Aggregate counts for duplicate gene names
tcga %>% filter(!grepl("^MT-", gene_name)) %>% filter(gene_type %in% good_types) %>% group_by(gene_name) %>% summarize_if(is.numeric,sum) %>% as.data.frame() -> dedup.ch.type.tcga #20328 genes, 425 samples
rownames(dedup.ch.type.tcga) <- dedup.ch.type.tcga$gene_name

# --- Filter and summarize Cell lines counts ---
# 1. Remove mitochondrial genes (MT-), 2. Filter by good gene type, 3. Aggregate counts for duplicate gene names
cell_lines %>% filter(!grepl("^MT-", gene_name)) %>% filter(gene_name %in% rownames(dedup.ch.type.tcga)) %>% group_by(gene_name) %>% summarize_if(is.numeric,sum) %>% as.data.frame() -> dedup.ch.type.cell_lines #20328 genes, 43 samples
rownames(dedup.ch.type.cell_lines) <- dedup.ch.type.cell_lines$gene_name

# Unite counts
total_counts <- as.data.frame(cbind(as.matrix(dedup.ch.type.tcga), as.matrix(dedup.ch.type.cell_lines))) #20328 genes, 468 samples
rownames(total_counts) <- total_counts$gene_name
all_counts$gene_name <- NULL

### PART 3. Metadata preparation and synchronization
# Prepare TCGA metadata
tumor_clusters <- tumor_clusters %>% select(Subtype, name2) %>% rename(cluster = Subtype, sample_id = name2) %>% mutate(protocol = "poly-A") 

# Prepare Cell lines metadata
cell_lines_clusters <- cell_lines_clusters %>% tibble::rownames_to_column("sample_id") %>% rename(cluster = Cell.line) 

# Merge metadata
total_meta <- merge(cell_lines_clusters, tumor_clusters, all = TRUE)
total_meta <- total_meta[!is.na(total_meta$protocol) & total_meta$protocol != "", ]

# Synchronize samples (counts and meta)
setdiff(total_meta$sample_id, colnames(total_counts))
total_meta$sample_id <- trimws(total_meta$sample_id)
colnames(total_counts) <- trimws(colnames(total_counts))
total_counts <- total_counts[, colnames(total_counts) %in% total_meta$sample_id]

# Removing samples that were visually identified as wrong outliers on PCA graphs created in earlier iterations of the analysis
samples_to_remove <- c(
  "SRR8060845",
  "SRR8060847",
  "SRR8060850",
  "SRR8060854",
  "SRR9995061"
)

total_counts <- total_counts[, !colnames(total_counts) %in% samples_to_remove]
total_meta <- total_meta %>% filter(!sample_id %in% samples_to_remove)

# Correct specific sample protocol
total_meta <- total_meta %>% mutate(protocol = if_else(sample_id == "SRR11539513", "Total RNA", protocol))


## PART 4. Gene filtering and protocol bias ratio
# Ensure all count data is numeric for calculation
total_counts_numeric <- total_counts %>% mutate(across(everything(), as.numeric))

# --- Filter non-variable genes (variance > 0) ---
gene_variance <- apply(total_counts_numeric, 1, var)
genes_to_keep <- which(gene_variance > 0)
filtered_counts_numeric <- total_counts_numeric[genes_to_keep, ] #20242

# --- Calculate protocol bias ratio ---
# Used to identify genes with extreme differences in expression between sequencing protocols
polyA_sum <- rowSums(filtered_counts_numeric[ , total_meta$sample_id[total_meta$protocol == 'poly-A']], na.rm = TRUE) # poly_A samples (TCGA + some cell-lines) are selected and summed numerically (424)
total_sum <- rowSums(filtered_counts_numeric[ , total_meta$sample_id[total_meta$protocol == 'Total RNA']], na.rm = TRUE) # Total samples (some cell-lines) are selected and summed numerically (33)

# Calculate log10 ratio (Total RNA / poly-A). Add a small constant to avoid division/log of zero
zeros = log10((total_sum + 10^-10)/polyA_sum)
sum((zeros == Inf)) #404

# Identify genes that are not expressed in poly-A (result in Inf ratio)
genes_to_use <- names(which(log_protocol_ratio != Inf))

# --- Initial DGEList creation and filtering ---
Total.DGE <- DGEList(filtered_counts_numeric[names(which(zeros != Inf)), total_meta$sample_id], 
                     genes = genes_to_use,
                     samples = total_meta$sample_id,
                     group = unlist(total_meta$cluster, use.names = F))

# Filter genes with zero total expression across all samples
zeroGenesI <- which(rowSums(Total.DGE$counts) == 0)
Total.DGE <- Total.DGE[setdiff(1:dim(Total.DGE)[1], zeroGenesI), ]

# Filter genes with low expression (mean log CPM <= 0)
smallCpmGenesI <- which(apply(cpm(Total.DGE, log = T), 1, mean) <= 0) #7207
Total.DGE <- Total.DGE[setdiff(1:dim( Total.DGE)[1], smallCpmGenesI), ]


### PART 5. Normalization and batch correction
# TMM normalization
Total.DGE <- calcNormFactors(Total.DGE,method = 'TMM')

# Calculate normalized log-CPM (before batch correction)
nonbatched <- log1p(cpm(Total.DGE))

# Batch correction
# Removes the variation attributed to the 'protocol' factor from the log-CPM expression
debatched = removeBatchEffect(log1p(cpm(Total.DGE)),batch = total_meta$protocol)

# Save normalized log-CPM matrix
write.table(debatched, file = "debatched_without_5.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)


### PART 6. Visualization (PCA) and Outlier Analysis
# --- Principal Component Analysis (PCA) ---
# PCA is run on the *transpose* of the expression matrix (samples as rows)
# PCA before batch correction
pca_before <- prcomp(t(nonbatched), scale. = TRUE)

# PCA after batch correction
pca_after <- prcomp(t(debatched), scale. = TRUE)

# Auxiliary function for PCA plotting
plot_pca <- function(pca, meta, title) {
  df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    protocol = meta$protocol,
    cluster = meta$cluster
  )
  ggplot(df, aes(x = PC1, y = PC2, color = protocol, shape = cluster)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_minimal()
}

# Generate PCA plots to assess batch effect removal
plot_pca(pca_before, total_meta, "PCA before batch correction")
plot_pca(pca_after, total_meta, "PCA after batch correction")


# --- Outlier detection based on PCA (example) ---
# This section attempts to identify samples with PC1 coordinates below a specific threshold and refine the original datasets
sample_scores <- as.data.frame(pca_before$x)
pc1_coordinates <- sample_scores$PC1
names(pc1_coordinates) <- rownames(sample_scores)

# Set the threshold
threshold <- -100 

# Identify samples with PC1 coordinate below the threshold
highly_outlying_samples <- pc1_coordinates[(pc1_coordinates) < threshold]

# Create final data frame of outliers with associated metadata
outlying_df <- data.frame(
  Sample_ID = names(highly_outlying_samples), # length 41
  PC1_Coordinate = highly_outlying_samples
)

meta_for_join <- total_meta %>%
  select(sample_id, cluster, protocol)

outlying_df_final <- outlying_df %>%
  left_join(meta_for_join, by = c("Sample_ID" = "sample_id"))

# View the identified outliers
print(outlying_df_final)
