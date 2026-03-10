# Load GEO DEGs
gse <- read.csv("results/tables/GSE19804_significant_DEGs.csv", row.names = 1)

# Load TCGA DEGs
tcga <- read.csv("results/tables/TCGA_LUAD_significant_DEGs.csv", row.names = 1)

# Extract gene symbols
gse_genes <- rownames(gse)
tcga_genes <- rownames(tcga)

# Find overlap
common_genes <- intersect(gse_genes, tcga_genes)

write.csv(common_genes, "results/tables/Common_DEGs_GSE_TCGA.csv", row.names = FALSE)

cat("Number of overlapping DEGs:", length(common_genes), "\n")
