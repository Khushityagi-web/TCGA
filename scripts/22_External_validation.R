# 22_External_validation.R
library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(AnnotationDbi)
library(ggplot2)

gse <- getGEO(filename="data_processed/GSE31210_series_matrix.txt.gz")
expr <- exprs(gse)
expr <- log2(expr + 1)  # FIX
pheno <- pData(gse)

group <- ifelse(grepl("tumor", pheno$characteristics_ch1), "Tumor", "Normal")
group <- factor(group)
cat("Groups:\n"); print(table(group))

design <- model.matrix(~ group)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, number=Inf, adjust.method="BH")

symbols <- mapIds(hgu133plus2.db, keys=rownames(results),
                  column="SYMBOL", keytype="PROBEID", multiVals="first")
results$GeneSymbol <- symbols
results <- results[!is.na(results$GeneSymbol), ]
results <- results[!duplicated(results$GeneSymbol), ]

sig <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
cat("GSE31210 significant DEGs:", nrow(sig), "\n")

write.csv(results, "results/tables/GSE31210_full_DEG.csv", row.names=FALSE)
write.csv(sig, "results/tables/GSE31210_significant_DEGs.csv", row.names=FALSE)

hub_genes <- c("DLGAP5","CCNA2","CCNB1","CDK1","KIF20A",
               "NCAPG","BUB1B","ASPM","KIF11","MAD2L1")

concordant <- read.csv("results/tables/direction_consistency.csv")
concordant <- concordant[concordant$Concordant==TRUE, ]

cat("\nHub gene direction consistency in GSE31210:\n")
validation_results <- list()
for (gene in hub_genes) {
  row <- results[results$GeneSymbol==gene, ]
  orig <- concordant[concordant$GeneSymbol==gene, "logFC_GSE"]
  if (nrow(row)>0 & length(orig)>0) {
    consistent <- sign(row$logFC)==sign(orig)
    sig_flag <- row$adj.P.Val < 0.05
    cat(gene, "- logFC:", round(row$logFC,3),
        "| adjP:", round(row$adj.P.Val,4),
        "| Consistent:", consistent,
        "| Significant:", sig_flag, "\n")
    validation_results[[gene]] <- data.frame(
      Gene=gene,
      logFC_GSE31210=round(row$logFC,3),
      adjP=round(row$adj.P.Val,4),
      Consistent=consistent,
      Significant=sig_flag)
  }
}

val_df <- do.call(rbind, validation_results)
write.csv(val_df, "results/tables/Hub_genes_GSE31210_validation.csv", row.names=FALSE)

overlap <- intersect(sig$GeneSymbol, concordant$GeneSymbol)
cat("\nOverlap with original concordant DEGs:", length(overlap), "\n")

results$Significant <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, "Yes", "No")
p <- ggplot(results, aes(x=logFC, y=-log10(adj.P.Val), color=Significant)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  labs(title="Volcano Plot - GSE31210 (External Validation)",
       x="Log2 Fold Change", y="-Log10 Adjusted P-value")
ggsave("results/figures/GSE31210_volcano.png", p, width=8, height=6)
cat("Done.\n")
