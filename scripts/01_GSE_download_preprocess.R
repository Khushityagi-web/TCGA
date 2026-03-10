library(GEOquery)
library(affy)
library(limma)

# Set working directory safely
setwd(getwd())

# Download GSE19804
gse <- getGEO("GSE19804", GSEMatrix = TRUE)
exprSet <- exprs(gse[[1]])

# Extract phenotype data
pheno <- pData(gse[[1]])

# Save raw expression
write.csv(exprSet, "data_processed/GSE19804_expression_raw.csv")
write.csv(pheno, "data_processed/GSE19804_pheno.csv")

cat("GSE19804 downloaded and saved\n")
