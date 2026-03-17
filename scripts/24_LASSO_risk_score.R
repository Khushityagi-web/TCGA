# 24_LASSO_risk_score.R
library(SummarizedExperiment)
library(glmnet)
library(survival)
library(survminer)
library(survcomp)
library(ggplot2)

se         <- readRDS("data_processed/TCGA_LUAD_expression.rds")
surv_df    <- read.csv("results/tables/Survival_data.csv")
counts     <- assay(se, "fpkm_unstrand")
gene_names <- rowData(se)$gene_name

hub_genes <- c("DLGAP5","CCNA2","CCNB1","CDK1","KIF20A",
               "NCAPG","BUB1B","ASPM","KIF11","MAD2L1")

# Build expression matrix
expr_mat <- data.frame(sample=colnames(counts))
for (gene in hub_genes) {
  idx <- which(gene_names==gene)
  if (length(idx)==0) next
  expr_mat[[gene]] <- log2(as.numeric(counts[idx[1],])+1)
}

df <- merge(surv_df, expr_mat, by="sample")
df <- df[!is.na(df$OS_time) & df$OS_time>0, ]
cat("Samples for LASSO:", nrow(df), "\n")

# LASSO Cox
x <- as.matrix(df[, hub_genes])
y <- Surv(df$OS_time, df$OS_status)
set.seed(42)
cv_fit <- cv.glmnet(x, y, family="cox", alpha=1, nfolds=10)

# Coefficients
coefs <- coef(cv_fit, s="lambda.min")
coef_df <- data.frame(
  Gene        = rownames(coefs),
  Coefficient = as.numeric(coefs)
)
coef_df <- coef_df[coef_df$Coefficient != 0, ]
cat("\nSelected genes by LASSO:\n")
print(coef_df)
write.csv(coef_df, "results/tables/LASSO_coefficients.csv", row.names=FALSE)

# Risk score
selected_genes <- coef_df$Gene
risk_score <- as.numeric(x[, selected_genes, drop=FALSE] %*% coef_df$Coefficient)
df$risk_score <- risk_score
df$risk_group <- ifelse(risk_score >= median(risk_score), "High", "Low")
cat("\nRisk group distribution:\n")
print(table(df$risk_group))

# KM plot
fit <- survfit(Surv(OS_time, OS_status) ~ risk_group, data=df)
p <- ggsurvplot(fit, data=df, pval=TRUE, risk.table=TRUE,
                title="LASSO Risk Score - Overall Survival",
                legend.labs=c("High Risk","Low Risk"),
                palette=c("#E7298A","#1B9E77"),
                xlab="Days", ylab="Survival Probability")
png("results/figures/LASSO_KM_risk_score.png", width=1000, height=800)
print(p, newpage=FALSE)
dev.off()
cat("KM plot saved\n")

# Risk score distribution
p2 <- ggplot(df, aes(x=reorder(sample, risk_score), y=risk_score, fill=risk_group)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("High"="#E7298A","Low"="#1B9E77")) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(title="Risk Score Distribution", x="Patients", y="Risk Score")
ggsave("results/figures/LASSO_risk_distribution.png", p2, width=10, height=5)
cat("Distribution plot saved\n")

# C-index
cindex <- concordance.index(
  x      = df$risk_score,
  surv.t = df$OS_time,
  surv.e = df$OS_status
)
cat("\nC-index:", round(cindex$c.index, 3), "\n")
cat("95% CI:", round(cindex$lower, 3), "-", round(cindex$upper, 3), "\n")

write.csv(df[, c("sample","risk_score","risk_group")],
          "results/tables/Patient_risk_scores.csv", row.names=FALSE)
cat("\nLASSO risk score analysis complete.\n")
