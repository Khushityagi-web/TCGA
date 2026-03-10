library(TCGAbiolinks)

query_maf <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query_maf)

maf <- GDCprepare(query_maf)

saveRDS(maf, "data_processed/TCGA_LUAD_maf.rds")

cat("TCGA mutation data downloaded\n")
