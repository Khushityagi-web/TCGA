# 14_Mutation_analysis.R
library(maftools)

maf_data <- readRDS("data_processed/TCGA_LUAD_maf.rds")

# Convert to maftools object
maf <- read.maf(maf = maf_data)

# Summary stats
write.csv(getSampleSummary(maf), "results/tables/MAF_sample_summary.csv", row.names=FALSE)
write.csv(getGeneSummary(maf),   "results/tables/MAF_gene_summary.csv",   row.names=FALSE)

# Plots
png("results/figures/MAF_summary.png", width=1200, height=800)
plotmafSummary(maf, rmOutlier=TRUE, addStat="median")
dev.off()

png("results/figures/Oncoplot_top20.png", width=1200, height=800)
oncoplot(maf, top=20)
dev.off()

# Get mutated genes list
gene_summary <- getGeneSummary(maf)
mutated_genes <- gene_summary$Hugo_Symbol

# Overlap with concordant DEGs
concordant <- read.csv("results/tables/direction_consistency.csv")
concordant <- concordant[concordant$Concordant == TRUE, ]

overlap_mut_deg <- intersect(mutated_genes, concordant$GeneSymbol)
write.csv(data.frame(GeneSymbol=overlap_mut_deg), 
          "results/tables/Mutated_Concordant_DEGs.csv", row.names=FALSE)

cat("Total mutated genes:", length(mutated_genes), "\n")
cat("Concordant DEGs also mutated:", length(overlap_mut_deg), "\n")
