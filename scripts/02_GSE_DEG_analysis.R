library(limma)

# Load data
expr <- read.csv("data_processed/GSE19804_expression_raw.csv", row.names = 1)
pheno <- read.csv("data_processed/GSE19804_pheno.csv", row.names = 1)

# Define groups using source_name_ch1
group <- ifelse(
  pheno$source_name_ch1 == "frozen tissue of primary tumor",
  "Tumor",
  "Normal"
)

group <- factor(group)

cat("Group distribution:\n")
print(table(group))

# Design matrix
design <- model.matrix(~ group)

# Fit model
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# Extract full results
results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

write.csv(results, "results/tables/GSE19804_full_DEG_table.csv")

# Filter significant DEGs
sig <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
write.csv(sig, "results/tables/GSE19804_significant_DEGs.csv")

cat("Total significant DEGs:", nrow(sig), "\n")
