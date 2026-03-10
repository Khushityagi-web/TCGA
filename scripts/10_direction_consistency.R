# Load symbol-mapped DEGs
gse <- read.csv("results/tables/GSE19804_DEGs_with_symbols.csv")
tcga <- read.csv("results/tables/TCGA_LUAD_DEGs_with_symbols.csv")

# Remove duplicates
gse <- gse[!duplicated(gse$GeneSymbol), ]
tcga <- tcga[!duplicated(tcga$GeneSymbol), ]

# Merge by gene symbol
merged <- merge(gse[, c("GeneSymbol", "logFC")],
                tcga[, c("GeneSymbol", "log2FoldChange")],
                by = "GeneSymbol")

colnames(merged) <- c("GeneSymbol", "logFC_GSE", "logFC_TCGA")

# Determine concordance
merged$Direction_GSE  <- ifelse(merged$logFC_GSE > 0, "Up", "Down")
merged$Direction_TCGA <- ifelse(merged$logFC_TCGA > 0, "Up", "Down")

merged$Concordant <- merged$Direction_GSE == merged$Direction_TCGA

cat("Total overlapping genes:", nrow(merged), "\n")
cat("Concordant genes:", sum(merged$Concordant), "\n")
cat("Discordant genes:", sum(!merged$Concordant), "\n")
