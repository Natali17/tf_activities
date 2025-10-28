### PART 1. Preprocessing
# Set the directory with the count_file
setwd('/Users/natal/PycharmProjects/contests/promoters/october2025')

# Load RNA expression counts for cell lines
counts <- read.delim("cell_lines_total_counts.tsv", row.names = 1, check.names = FALSE)

# Remove gene name column
counts$gene_name <- NULL

# Remove version numbers from Ensembl IDs and keep unique genes
counts$Ensembl_ID <- sub("\\..*", "", rownames(counts))
counts_filtered <- counts[!duplicated(counts$Ensembl_ID), ]
rownames(counts_filtered) <- counts_filtered$Ensembl_ID
counts_filtered$Ensembl_ID <- NULL

# Remove samples with unknown protocol
cols_to_remove <- c("SRR12144191", "SRR19909308", "SRR19909310", "SRR19909307")
counts_filtered <- counts_filtered[, !(colnames(counts_filtered) %in% cols_to_remove)]

# Load sample metadata
sample_info <- read.csv("Cell_lines_protocols.csv", header = TRUE)
sample_info <- sample_info[sample_info$protocol != "", ]
sample_info$SRR <- trimws(sample_info$SRR)

# Check that metadata names match the expression data
all(sample_info$SRR %in% colnames(counts_filtered))
#setdiff(colnames(counts_filtered), sample_info$SRR)

# Reorder metadata to match expression data columns
sample_info <- sample_info[match(colnames(counts_filtered), sample_info$SRR), ]

### PART 2. TMM-normalization and visualization
#BiocManager::install('edgeR')
library(edgeR)

# Define sample groups (factor levels) for analysis
HepD <- factor(sample_info$Cell.line, levels=c("Hep3B","HepG2", "HUH7"))

# Create DGEList object for downstream analysis
DGE <- DGEList(counts = counts_filtered, group = HepD)

# Filter out genes with very low expression
keep <- filterByExpr(DGE, group = HepD)
DGE <- DGE[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes using TMM
DGE <- calcNormFactors(DGE, method = "TMM")

# Convert to log2 CPM with prior count to avoid log(0)
log2_cpm <- cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

# Inspect samples and group counts
head(DGE$samples)
table(DGE$samples$group)

# Save normalized log2 CPM matrix
write.table(log2_cpm, file = "cell_lines_log2tmm.tsv",
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
                     CellLine = HepD)
ggplot(pca_df, aes(x = PC1, y = PC2, color = CellLine)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: различия между клеточными линиями (TMM-normalized)")

##########################################################
### PCA не очень кластеризует, сделаем фильтрацию и повторим визуализацию
# Выбираем, например, топ 5000 генов по дисперсии
var_genes <- apply(log2_cpm, 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:5000]
log2_cpm_top <- log2_cpm[top_genes, ]
pca_top <- prcomp(t(log2_cpm_top))

pca_df_top <- data.frame(PC1 = pca_top$x[,1],
                     PC2 = pca_top$x[,2],
                     CellLine = HepD)

ggplot(pca_df_top, aes(x = PC1, y = PC2, color = CellLine)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: различия между клеточными линиями (TMM-normalized)")


### 
install.packages('umap')
library(umap)
umap_res <- umap(t(log2_cpm))
plot(umap_res$layout, col = as.factor(DGE$samples$group), pch=19)
library(pheatmap)
pheatmap(log2_cpm_top, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, annotation_col = data.frame(CellLine=DGE$samples$group))

annotation_col <- data.frame(CellLine = DGE$samples$group)
rownames(annotation_col) <- colnames(log2_cpm_top)

pheatmap(log2_cpm_top,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         scale = "row")  # стандартизируем по рядам


# 4️⃣ Проверить batch-эффекты

#Если данные получены из разных экспериментов или протоколов, они могут затмевать биологические различия. Тогда нужно:
  
 # включить batch в дизайн при дифф. экспрессии (model.matrix(~ batch + CellLine)),

#или корректировать данные с помощью limma::removeBatchEffect() для PCA/визуализации.
# Если есть batch (например, protocol или dataset), можно посмотреть влияние:
if("protocol" %in% colnames(sample_info)){
  batch <- factor(sample_info$protocol)
  # Коррекция для визуализации
  library(limma)
  log2_cpm_batch_corrected <- removeBatchEffect(log2_cpm_top, batch = batch)
  
  pheatmap(log2_cpm_batch_corrected,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           annotation_col = annotation_col,
           scale = "row",
           main = "Batch-corrected Heatmap")
}
