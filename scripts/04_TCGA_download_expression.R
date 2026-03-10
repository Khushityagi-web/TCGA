library(TCGAbiolinks)

setwd(getwd())

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)

data <- GDCprepare(query)

saveRDS(data, "data_processed/TCGA_LUAD_expression.rds")

cat("TCGA RNA-seq data downloaded and saved\n")
