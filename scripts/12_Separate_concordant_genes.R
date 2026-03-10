# 12_Separate_concordant_genes.R
gse <- read.csv("results/tables/GSE19804_DEGs_with_symbols.csv")
tcga <- read.csv("results/tables/TCGA_LUAD_DEGs_with_symbols.csv")

gse <- gse[!duplicated(gse$GeneSymbol), ]
tcga <- tcga[!duplicated(tcga$GeneSymbol), ]

merged <- merge(gse[, c("GeneSymbol", "logFC")],
                tcga[, c("GeneSymbol", "log2FoldChange")],
                by = "GeneSymbol")
colnames(merged) <- c("GeneSymbol", "logFC_GSE", "logFC_TCGA")

merged$Direction_GSE  <- ifelse(merged$logFC_GSE > 0, "Up", "Down")
merged$Direction_TCGA <- ifelse(merged$logFC_TCGA > 0, "Up", "Down")
merged$Concordant <- merged$Direction_GSE == merged$Direction_TCGA

# Save merged table this time
write.csv(merged, "results/tables/direction_consistency.csv", row.names=FALSE)

concordant <- merged[merged$Concordant == TRUE, ]
up <- concordant[concordant$Direction_GSE == "Up", "GeneSymbol"]
down <- concordant[concordant$Direction_GSE == "Down", "GeneSymbol"]

write.csv(data.frame(GeneSymbol=up), "results/tables/Concordant_UP_genes.csv", row.names=FALSE)
write.csv(data.frame(GeneSymbol=down), "results/tables/Concordant_DOWN_genes.csv", row.names=FALSE)

cat("Concordant UP:", length(up), "\n")
cat("Concordant DOWN:", length(down), "\n")
