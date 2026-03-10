# 18_KM_top_survival_genes.R
library(SummarizedExperiment)
library(survival)
library(survminer)
library(ggplot2)

se         <- readRDS("data_processed/TCGA_LUAD_expression.rds")
surv_df    <- read.csv("results/tables/Survival_data.csv")
counts     <- assay(se, "fpkm_unstrand")
gene_names <- rowData(se)$gene_name

res   <- read.csv("results/tables/Survival_screen_all.csv")
top10 <- head(res[res$FDR < 0.05, "Gene"], 10)

for (gene in top10) {
  idx <- which(gene_names == gene)
  if (length(idx) == 0) next

  expr_vals        <- log2(as.numeric(counts[idx[1], ]) + 1)
  names(expr_vals) <- colnames(counts)

  df       <- surv_df
  df$expr  <- expr_vals[df$sample]
  df$group <- ifelse(df$expr >= median(df$expr, na.rm=TRUE), "High", "Low")

  fit <- survfit(Surv(OS_time, OS_status) ~ group, data=df)
  hr  <- res[res$Gene == gene, "HR"]

  p <- ggsurvplot(fit, data=df, pval=TRUE, risk.table=TRUE,
                  title=paste0(gene, " (HR=", hr, ")"),
                  legend.labs=c("High","Low"),
                  palette=c("#E7298A","#1B9E77"),
                  xlab="Days", ylab="Survival Probability")

  # Correct way to save ggsurvplot
  png(paste0("results/figures/KM_", gene, ".png"), width=1000, height=800)
  print(p)
  dev.off()
  cat("Saved:", gene, "\n")
}
