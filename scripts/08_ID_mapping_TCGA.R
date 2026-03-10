library(org.Hs.eg.db)
library(AnnotationDbi)

# Load TCGA DEGs
tcga <- read.csv("results/tables/TCGA_LUAD_significant_DEGs.csv", row.names=1)

# Remove version numbers from Ensembl IDs
ensembl_ids <- gsub("\\..*", "", rownames(tcga))

# Map Ensembl IDs to Gene Symbols
symbols <- mapIds(org.Hs.eg.db,
                  keys = ensembl_ids,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

tcga$GeneSymbol <- symbols

# Remove unmapped genes
tcga_clean <- tcga[!is.na(tcga$GeneSymbol), ]

write.csv(tcga_clean, "results/tables/TCGA_LUAD_DEGs_with_symbols.csv")

cat("Mapped TCGA genes:", nrow(tcga_clean), "\n")
