library(DESeq2)
library(SummarizedExperiment)

# Load data
se <- readRDS("data_processed/TCGA_LUAD_expression.rds")

# Extract counts
counts <- assay(se, "unstranded")

# Extract metadata
meta <- colData(se)

# Identify tumor vs normal
group <- ifelse(substr(meta$shortLetterCode, 1, 2) == "NT", "Normal", "Tumor")
group <- factor(group)

cat("Group distribution:\n")
print(table(group))

# Filter low count genes
keep <- rowSums(counts >= 10) >= 10
counts <- counts[keep, ]

# Create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(group),
  design = ~ group
)

dds <- DESeq(dds)

res <- results(dds)

resOrdered <- res[order(res$padj), ]

# Save full results
write.csv(as.data.frame(resOrdered), "results/tables/TCGA_LUAD_full_DEG_table.csv")

# Filter significant DEGs
sig <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig), "results/tables/TCGA_LUAD_significant_DEGs.csv")

cat("Total significant TCGA DEGs:", nrow(sig), "\n")
