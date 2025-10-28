### PART 1. Preprocessing
# Set the directory with the count_file
setwd('/Users/natal/PycharmProjects/contests/promoters/october2025')

# Load tumor cluster metadata
tumor_clusters <- read.csv("TCGA_only_labeled_meta.csv", header = TRUE, row.names = 1)

# Replace dots with hyphens in sample names for consistency
tumor_clusters$name <- gsub("\\.", "-", tumor_clusters$name)

# Load TCGA gene expression counts
tcga_total_counts <- read.table('tcga_total_counts.tsv', header = TRUE, row.names = 1)

# Remove unnecessary annotation columns
cols_to_remove <- c('gene_name', 'gene_type')
tcga_total_counts <- tcga_total_counts[ , !(colnames(tcga_total_counts) %in% cols_to_remove)]

# Remove version numbers from Ensembl IDs and keep unique genes
tcga_total_counts$Ensembl_ID <- sub("\\..*", "", rownames(tcga_total_counts))
tcga_total_counts <- tcga_total_counts[!duplicated(tcga_total_counts$Ensembl_ID), ]
rownames(tcga_total_counts) <- tcga_total_counts$Ensembl_ID
tcga_total_counts$Ensembl_ID <- NULL
# Check that metadata names match the expression data
all(tumor_clusters$name2 %in% colnames(tcga_total_counts))

# Reorder metadata to match expression data columns
tumor_clusters <- tumor_clusters[match(colnames(tcga_total_counts), tumor_clusters$name2), ]

### PART 2. TMM-normalization and visualization
#BiocManager::install('edgeR')
library(edgeR)

# Define sample groups (factor levels) for analysis
states <- factor(tumor_clusters$Subtype, levels=c("DDF Tumor","NAT", "Standard Tumor"))

# Create DGEList object for downstream analysis
DGE <- DGEList(counts = tcga_total_counts, group = states)

# Filter out genes with very low expression
keep <- filterByExpr(DGE, group = states)
DGE <- DGE[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes using TMM
DGE <- calcNormFactors(DGE, method = "TMM")

# Convert to log2 CPM with prior count to avoid log(0)
log2_cpm <- cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

# Inspect samples and group counts
head(DGE$samples)
table(DGE$samples$group)

# Save normalized log2 CPM matrix
write.table(log2_cpm, file = "tcga_log2tmm.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

### PART 3. ADDITIONAL INSPECTION
library(ggplot2)

# Perform PCA on transposed log2 CPM matrix
pca <- prcomp(t(log2_cpm))
pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     States = states)
ggplot(pca_df, aes(x = PC1, y = PC2, color = States)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: различия между DDF, NAT и Standart Tumor")
###