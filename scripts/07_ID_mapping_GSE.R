library(hgu133plus2.db)
library(AnnotationDbi)

# Load GSE DEGs
gse <- read.csv("results/tables/GSE19804_significant_DEGs.csv", row.names=1)

probe_ids <- rownames(gse)

# Map probe IDs to gene symbols
symbols <- mapIds(hgu133plus2.db,
                  keys = probe_ids,
                  column = "SYMBOL",
                  keytype = "PROBEID",
                  multiVals = "first")

gse$GeneSymbol <- symbols

# Remove NA mappings
gse_clean <- gse[!is.na(gse$GeneSymbol), ]

write.csv(gse_clean, "results/tables/GSE19804_DEGs_with_symbols.csv")

cat("Mapped GSE genes:", nrow(gse_clean), "\n")
