# 16_Gene_survival.R
library(SummarizedExperiment)
library(survival)
library(survminer)
library(ggplot2)

se      <- readRDS("data_processed/TCGA_LUAD_expression.rds")
surv_df <- read.csv("results/tables/Survival_data.csv")

counts     <- assay(se, "fpkm_unstrand")
gene_names <- rowData(se)$gene_name

concordant <- read.csv("results/tables/direction_consistency.csv")
concordant <- concordant[concordant$Concordant == TRUE, ]
concordant$meanFC <- (abs(concordant$logFC_GSE) + abs(concordant$logFC_TCGA)) / 2
mutated  <- read.csv("results/tables/Mutated_Concordant_DEGs.csv")
priority <- concordant[concordant$GeneSymbol %in% mutated$GeneSymbol, ]
top10    <- head(priority[order(-priority$meanFC), "GeneSymbol"], 10)
cat("Top 10 priority genes:", paste(top10, collapse=", "), "\n")

results_list <- list()

for (gene in top10) {
  idx <- which(gene_names == gene)
  if (length(idx) == 0) { cat("Not found:", gene, "\n"); next }

  # Log2 transform - this is the fix
  expr_vals <- log2(as.numeric(counts[idx[1], ]) + 1)
  names(expr_vals) <- colnames(counts)

  common <- intersect(surv_df$sample, names(expr_vals))
  if (length(common) < 50) next

  df       <- surv_df[surv_df$sample %in% common, ]
  df$expr  <- expr_vals[df$sample]
  df$group <- ifelse(df$expr >= median(df$expr, na.rm=TRUE), "High", "Low")

  cox  <- coxph(Surv(OS_time, OS_status) ~ expr, data=df)
  pval <- summary(cox)$coefficients[5]
  hr   <- summary(cox)$coefficients[2]

  results_list[[gene]] <- data.frame(Gene=gene, HR=round(hr,3), P_value=round(pval,4))

  # Plot regardless of significance for top genes
  fit <- survfit(Surv(OS_time, OS_status) ~ group, data=df)
  p   <- ggsurvplot(fit, data=df, pval=TRUE, risk.table=TRUE,
                    title=paste(gene, "- OS Stratification (High vs Low)"),
                    legend.labs=c("High","Low"),
                    palette=c("#E7298A","#1B9E77"))
  ggsave(paste0("results/figures/KM_", gene, ".png"), print(p), width=10, height=7)
}

survival_results <- do.call(rbind, results_list)
write.csv(survival_results, "results/tables/Gene_survival_results.csv", row.names=FALSE)
cat("\nSurvival results:\n")
print(survival_results)
cat("\nGenes significant (p<0.05):", sum(survival_results$P_value < 0.05, na.rm=TRUE), "\n")
