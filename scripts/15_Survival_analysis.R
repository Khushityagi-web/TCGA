# 15_Survival_analysis.R
library(SummarizedExperiment)
library(survival)
library(survminer)
library(ggplot2)

se <- readRDS("data_processed/TCGA_LUAD_expression.rds")
clinical <- colData(se)

surv_df <- data.frame(
  sample     = clinical$barcode,
  days       = as.numeric(clinical$days_to_death),
  days_fu    = as.numeric(clinical$days_to_last_follow_up),
  vital      = clinical$vital_status,
  stringsAsFactors = FALSE
)

surv_df$OS_time   <- ifelse(!is.na(surv_df$days), surv_df$days, surv_df$days_fu)
surv_df$OS_status <- ifelse(surv_df$vital == "Dead", 1, 0)
surv_df <- surv_df[!is.na(surv_df$OS_time) & surv_df$OS_time > 0, ]

# Overall KM curve
fit <- survfit(Surv(OS_time, OS_status) ~ 1, data=surv_df)
p <- ggsurvplot(fit, data=surv_df, risk.table=TRUE,
                title="Overall Survival - TCGA LUAD",
                xlab="Days", ylab="Survival Probability")
ggsave("results/figures/KM_overall.png",
       print(p), width=10, height=7)

# Top 5 concordant genes by mean |logFC|
concordant <- read.csv("results/tables/direction_consistency.csv")
concordant <- concordant[concordant$Concordant==TRUE, ]
concordant$meanFC <- (abs(concordant$logFC_GSE) + abs(concordant$logFC_TCGA)) / 2
top5 <- head(concordant[order(-concordant$meanFC), "GeneSymbol"], 5)
cat("Top 5 genes for survival analysis:", paste(top5, collapse=", "), "\n")

write.csv(surv_df, "results/tables/Survival_data.csv", row.names=FALSE)
cat("Survival data saved:", nrow(surv_df), "patients\n")
