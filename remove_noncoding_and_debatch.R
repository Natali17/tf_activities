library(dplyr)
library(limma)
library(edgeR)

### PART 1. Preprocessing
# Set the directory with the necessary files
setwd("//wsl.localhost/Ubuntu/home/natal/maradoner_project")

# Load TCGA and cell lines expression data
tcga <- read.table("C:/users/natal/PycharmProjects/contests/promoters/october2025/tcga_total_counts.tsv", header = 1, row.names = 1)
cell_lines <- read.table("C:/users/natal/PycharmProjects/contests/promoters/october2025/cell_lines_total_counts.tsv", header = 1, row.names = 1)
  
# Define gene types that should be remained
good_types <- c('IG_C_gene','IG_V_gene','protein_coding','TR_V_gene','transcribed_unitary_pseudogene')

# Filter expression data
tcga %>% filter(gene_type %in% good_types) %>% group_by(gene_name) %>% summarize_if(is.numeric,sum) %>% as.data.frame() -> dedup.ch.type.tcga
rownames(dedup.ch.type.tcga) <- dedup.ch.type.tcga$gene_name

cell_lines %>% filter(gene_name %in% rownames(dedup.ch.type.tcga)) %>% group_by(gene_name) %>% summarize_if(is.numeric,sum) %>% as.data.frame() -> dedup.ch.type.cell_lines
rownames(dedup.ch.type.cell_lines) <- dedup.ch.type.cell_lines$gene_name

# Unite tcga data with cell_line data
all_counts <- merge(dedup.ch.type.tcga, dedup.ch.type.cell_lines, by.x = "gene_name", by.y = "gene_name")
rownames(all_counts) <- all_counts$gene_name
all_counts$gene_name <- NULL

# Take only the genes with a total expression sum > 0
# expressed_genes <- rownames(dedup.ch.type.tcga)[rowSums(dedup.ch.type.tcga[ , 2:ncol(dedup.ch.type.tcga)], na.rm = TRUE) > 0]

# Load expression metadata
tumor_clusters <- read.csv("C:/Users/natal/PycharmProjects/contests/promoters/october2025/TCGA_only_labeled_meta.csv", header = TRUE, row.names = 1)
cell_lines_clusters <- read.csv("C:/Users/natal/PycharmProjects/contests/promoters/october2025/Cell_lines_protocols.csv", header = TRUE, row.names = 1)

# 1. Суммируем по cell_lines (только poly-A колонки)
cell_lines_numeric <- dedup.ch.type.cell_lines[, sapply(dedup.ch.type.cell_lines, is.numeric)]
cell_lines_sum <- as.data.frame(rowSums(cell_lines_numeric[, rownames(cell_lines_clusters)[cell_lines_clusters$protocol == 'poly-A']], na.rm = TRUE
))
cell_lines_sum$expression <- cell_lines_sum$`rowSums(cell_lines_numeric[, rownames(cell_lines_clusters)[cell_lines_clusters$protocol == "poly-A"]], na.rm = TRUE)`
cell_lines_sum$`rowSums(cell_lines_numeric[, rownames(cell_lines_clusters)[cell_lines_clusters$protocol == "poly-A"]], na.rm = TRUE)` <- NULL
# 2. Суммируем по TCGA (все колонки)
dedup.ch.type.tcga[, 2:ncol(dedup.ch.type.tcga)] <- lapply(dedup.ch.type.tcga[, 2:ncol(dedup.ch.type.tcga)], function(x) as.numeric(as.character(x)))
tcga_numeric <- dedup.ch.type.tcga[, sapply(dedup.ch.type.tcga, is.numeric)]
tcga_sum <- as.data.frame(rowSums(tcga_numeric, na.rm = TRUE))
tcga_sum$expression <- tcga_sum$`rowSums(tcga_numeric, na.rm = TRUE)`
tcga_sum$`rowSums(tcga_numeric, na.rm = TRUE)` <- NULL

# 3. Объединяем суммы по генам
# Сначала согласуем порядок генов
common_genes <- intersect(names(cell_lines_sum), names(tcga_sum))
zeros_AAA <- cell_lines_sum[common_genes] + tcga_sum[common_genes]

# 4. Суммируем по cell_lines (только total RNA колонки)
setdiff(
  rownames(cell_lines_clusters)[cell_lines_clusters$protocol == "Total RNA"],
  colnames(cell_lines_numeric)
)
rownames(cell_lines_clusters) <- trimws(rownames(cell_lines_clusters))

zeros_total <- as.data.frame(rowSums(dedup.ch.type.cell_lines[,rownames(cell_lines_clusters)[cell_lines_clusters$protocol == 'Total RNA']])) # Выбираются сэмплы total (some cell-lines) и суммируются погенно
zeros_total$expression <- zeros_total$`rowSums(dedup.ch.type.cell_lines[, rownames(cell_lines_clusters)[cell_lines_clusters$protocol == "Total RNA"]])`
zeros_total$`rowSums(dedup.ch.type.cell_lines[, rownames(cell_lines_clusters)[cell_lines_clusters$protocol == "Total RNA"]])` <- NULL

zeros = log10((zeros_total + 10^-10)/zeros_AAA) # для того чтобы убрать гены которые не экспрессируются в poly-A


tumor_clusters <- tumor_clusters %>% select(Subtype, name2) %>% rename(cluster = Subtype, sample_id = name2) %>% mutate(protocol = "poly-A") 
cell_lines_clusters <- cell_lines_clusters %>% tibble::rownames_to_column("sample_id") %>% rename(cluster = Cell.line) 

total_meta <- merge(cell_lines_clusters, tumor_clusters, all = TRUE)
total_meta <- total_meta[!is.na(total_meta$protocol) & total_meta$protocol != "", ]

setdiff(total_meta$sample_id, colnames(all_counts))
total_meta$sample_id <- trimws(total_meta$sample_id)
all_counts <- all_counts[, colnames(all_counts) %in% total_meta$sample_id]

### Part 2.TMM-normalization and removing batch effect
# Define sample groups (factor levels) for analysis
groups <- factor(total_meta$sample_id, levels=c("DDF Tumor","NAT", "Standard Tumor", "Hep3B", "HepG2", "HUH7"))

# Create DGEList object for downstream analysis
valid_genes <- rownames(zeros)[which(zeros$expression != Inf)]
valid_samples <- intersect(total_meta$sample_id, colnames(all_counts))


Total.DGE <- DGEList(counts  = all_counts[valid_genes, valid_samples], 
                     genes = valid_genes,
                     samples = valid_samples,
                     group = groups)
Total.DGE <- calcNormFactors(Total.DGE,method = 'TMM')

# Convert to log2 CPM with prior count to avoid log(0)
log_tmm <- cpm(Total.DGE, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1) 
debatched = removeBatchEffect(log_tmm,batch = total_meta$protocol)



# Save normalized log2 TMM matrix
write.table(debatched, file = "debatched_clusters.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

library(readr)
write_csv(total_meta, "total_meta.csv")


#
library(ggplot2)
sum(is.infinite(log_tmm)) 

# PCA ДО удаления batch-эффекта
pca_before <- prcomp(t(log_tmm), scale. = TRUE)

# PCA ПОСЛЕ удаления batch-эффекта
pca_after <- prcomp(t(debatched), scale. = TRUE)

# Вспомогательная функция для построения PCA-графика
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

# PCA-графики
plot_pca(pca_before, total_meta, "PCA before batch correction")
plot_pca(pca_after, total_meta, "PCA after batch correction")




