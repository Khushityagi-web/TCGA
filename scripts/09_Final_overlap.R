# Load cleaned symbol-based DEGs
gse <- read.csv("results/tables/GSE19804_DEGs_with_symbols.csv")
tcga <- read.csv("results/tables/TCGA_LUAD_DEGs_with_symbols.csv")

# Remove duplicates
gse_genes <- unique(gse$GeneSymbol)
tcga_genes <- unique(tcga$GeneSymbol)

# Find overlap
common <- intersect(gse_genes, tcga_genes)

write.csv(common, "results/tables/Common_DEGs_symbol_overlap.csv", row.names = FALSE)

cat("True overlapping genes:", length(common), "\n")
