# 20_Multivariate_Cox.R
library(SummarizedExperiment)
library(survival)

se         <- readRDS("data_processed/TCGA_LUAD_expression.rds")
surv_df    <- read.csv("results/tables/Survival_data.csv")
counts     <- assay(se, "fpkm_unstrand")
gene_names <- rowData(se)$gene_name

# Top 5 hub genes
hub_genes <- c("CCNB1", "CDK1", "BUB1B", "MAD2L1", "KIF11")

# Build expression matrix
expr_mat <- data.frame(sample=colnames(counts))
for (gene in hub_genes) {
  idx <- which(gene_names == gene)
  if (length(idx) == 0) next
  expr_mat[[gene]] <- log2(as.numeric(counts[idx[1], ]) + 1)
}

# Merge with survival
df <- merge(surv_df, expr_mat, by="sample")
df <- df[!is.na(df$OS_time) & df$OS_time > 0, ]
cat("Samples for multivariate analysis:", nrow(df), "\n")

# Multivariate Cox
formula_str <- paste("Surv(OS_time, OS_status) ~", paste(hub_genes, collapse=" + "))
cox_multi <- coxph(as.formula(formula_str), data=df)
summary_cox <- summary(cox_multi)

# Save results
coef_df <- as.data.frame(summary_cox$coefficients)
coef_df$Gene <- rownames(coef_df)
coef_df <- coef_df[, c("Gene", "exp(coef)", "Pr(>|z|)")]
colnames(coef_df) <- c("Gene", "HR", "P_value")
coef_df$FDR <- p.adjust(coef_df$P_value, method="BH")

write.csv(coef_df, "results/tables/Multivariate_Cox.csv", row.names=FALSE)
cat("\nMultivariate Cox results:\n")
print(coef_df)
cat("\nConcordance index:", summary_cox$concordance[1], "\n")
