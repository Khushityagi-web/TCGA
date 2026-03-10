# 21_Univariate_Cox_hubs.R
library(SummarizedExperiment)
library(survival)
library(ggplot2)

se         <- readRDS("data_processed/TCGA_LUAD_expression.rds")
surv_df    <- read.csv("results/tables/Survival_data.csv")
counts     <- assay(se, "fpkm_unstrand")
gene_names <- rowData(se)$gene_name

hub_genes <- c("NCAPG","CCNB1","CCNA2","BUB1B","ASPM",
               "KIF20A","CDK1","DLGAP5","KIF11","MAD2L1")

results <- list()
for (gene in hub_genes) {
  idx <- which(gene_names == gene)
  if (length(idx) == 0) next
  expr_vals        <- log2(as.numeric(counts[idx[1],]) + 1)
  names(expr_vals) <- colnames(counts)
  df               <- surv_df
  df$expr          <- expr_vals[df$sample]
  cox  <- coxph(Surv(OS_time, OS_status) ~ expr, data=df)
  s    <- summary(cox)
  results[[gene]] <- data.frame(
    Gene        = gene,
    HR          = round(s$coefficients[2], 3),
    CI_lower    = round(s$conf.int[3], 3),
    CI_upper    = round(s$conf.int[4], 3),
    P_value     = round(s$coefficients[5], 5)
  )
}

res_df <- do.call(rbind, results)
res_df$FDR <- p.adjust(res_df$P_value, method="BH")
res_df <- res_df[order(res_df$P_value), ]
write.csv(res_df, "results/tables/Univariate_Cox_hubs.csv", row.names=FALSE)

# Forest plot
res_df$Gene <- factor(res_df$Gene, levels=rev(res_df$Gene))
p <- ggplot(res_df, aes(x=HR, y=Gene, xmin=CI_lower, xmax=CI_upper)) +
  geom_pointrange(color="steelblue", size=0.8) +
  geom_vline(xintercept=1, linetype="dashed", color="red") +
  theme_minimal() +
  labs(title="Univariate Cox - Hub Genes",
       x="Hazard Ratio (95% CI)", y="Gene")
ggsave("results/figures/Forest_plot_hubs.png", p, width=8, height=6)

cat("\nUnivariate Cox results:\n")
print(res_df)
