# 23_Multivariate_Cox_clinical.R
library(SummarizedExperiment)
library(survival)

se         <- readRDS("data_processed/TCGA_LUAD_expression.rds")
surv_df    <- read.csv("results/tables/Survival_data.csv")
counts     <- assay(se, "fpkm_unstrand")
gene_names <- rowData(se)$gene_name
clinical   <- colData(se)

clin_df <- data.frame(
  sample = clinical$barcode,
  age    = as.numeric(clinical$age_at_index),
  gender = clinical$gender,
  stage  = clinical$ajcc_pathologic_stage,
  stringsAsFactors = FALSE
)

cat("Stage distribution:\n"); print(table(clin_df$stage))
cat("Gender distribution:\n"); print(table(clin_df$gender))

hub_genes <- c("DLGAP5","CCNA2","CDK1")
expr_mat  <- data.frame(sample=colnames(counts))
for (gene in hub_genes) {
  idx <- which(gene_names==gene)
  if (length(idx)==0) next
  expr_mat[[gene]] <- log2(as.numeric(counts[idx[1],])+1)
}

df <- merge(surv_df, expr_mat, by="sample")
df <- merge(df, clin_df, by="sample")
df <- df[!is.na(df$OS_time) & df$OS_time>0, ]

df$stage_simple <- ifelse(grepl("IV", df$stage), "IV",
                   ifelse(grepl("III", df$stage), "III",
                   ifelse(grepl("II", df$stage), "II", "I")))
df$stage_simple <- factor(df$stage_simple, levels=c("I","II","III","IV"))

cat("\nSamples for multivariate analysis:", nrow(df), "\n")
cat("Stage simplified:\n"); print(table(df$stage_simple))

cox <- coxph(Surv(OS_time, OS_status) ~
             DLGAP5 + CCNA2 + CDK1 + age + gender + stage_simple,
             data=df)
s <- summary(cox)
cat("\nMultivariate Cox results:\n")
print(round(s$coefficients[, c(1,2,5)], 4))
cat("\nConcordance:", round(s$concordance[1], 3), "\n")

write.csv(as.data.frame(s$coefficients),
          "results/tables/Multivariate_Cox_clinical.csv")
cat("Done.\n")
