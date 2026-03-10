# 17_Survival_screen.R
library(SummarizedExperiment)
library(survival)

se      <- readRDS("data_processed/TCGA_LUAD_expression.rds")
surv_df <- read.csv("results/tables/Survival_data.csv")
counts     <- assay(se, "fpkm_unstrand")
gene_names <- rowData(se)$gene_name

mutated    <- read.csv("results/tables/Mutated_Concordant_DEGs.csv")
concordant <- read.csv("results/tables/direction_consistency.csv")
concordant <- concordant[concordant$Concordant == TRUE, ]
genes_to_test <- intersect(mutated$GeneSymbol, concordant$GeneSymbol)
cat("Testing", length(genes_to_test), "genes\n")

results <- list()
for (gene in genes_to_test) {
  idx <- which(gene_names == gene)
  if (length(idx) == 0) next
  
  expr_vals        <- log2(as.numeric(counts[idx[1], ]) + 1)
  names(expr_vals) <- colnames(counts)
  
  df       <- surv_df
  df$expr  <- expr_vals[df$sample]
  if (sum(!is.na(df$expr)) < 50) next
  
  tryCatch({
    cox  <- coxph(Surv(OS_time, OS_status) ~ expr, data=df)
    pval <- summary(cox)$coefficients[5]
    hr   <- summary(cox)$coefficients[2]
    results[[gene]] <- data.frame(Gene=gene, HR=round(hr,3), P_value=round(pval,5))
  }, error=function(e) NULL)
}

res_df <- do.call(rbind, results)
res_df <- res_df[order(res_df$P_value), ]
res_df$FDR <- p.adjust(res_df$P_value, method="BH")

write.csv(res_df, "results/tables/Survival_screen_all.csv", row.names=FALSE)
cat("Significant (p<0.05):", sum(res_df$P_value < 0.05, na.rm=TRUE), "\n")
cat("Significant (FDR<0.05):", sum(res_df$FDR < 0.05, na.rm=TRUE), "\n")
cat("\nTop 20 survival genes:\n")
print(head(res_df, 20))
