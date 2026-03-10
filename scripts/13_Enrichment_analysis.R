# 13_Enrichment_analysis.R
library(enrichR)
library(ggplot2)

up <- read.csv("results/tables/Concordant_UP_genes.csv")$GeneSymbol
down <- read.csv("results/tables/Concordant_DOWN_genes.csv")$GeneSymbol

dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")

enrich_up   <- enrichr(up,   dbs)
enrich_down <- enrichr(down, dbs)

# Save tables
write.csv(enrich_up$GO_Biological_Process_2023,   "results/tables/GO_UP.csv",   row.names=FALSE)
write.csv(enrich_down$GO_Biological_Process_2023, "results/tables/GO_DOWN.csv", row.names=FALSE)
write.csv(enrich_up$KEGG_2021_Human,   "results/tables/KEGG_UP.csv",   row.names=FALSE)
write.csv(enrich_down$KEGG_2021_Human, "results/tables/KEGG_DOWN.csv", row.names=FALSE)

# Plot function
plot_enrichment <- function(df, title, filename) {
  df <- df[df$Adjusted.P.value < 0.05, ]
  df <- head(df[order(df$Adjusted.P.value), ], 15)
  df$Term <- factor(df$Term, levels=rev(df$Term))
  p <- ggplot(df, aes(x=-log10(Adjusted.P.value), y=Term)) +
    geom_bar(stat="identity", fill="steelblue") +
    theme_minimal() +
    labs(title=title, x="-log10(Adj.P)", y="")
  ggsave(filename, p, width=10, height=7)
}

plot_enrichment(enrich_up$GO_Biological_Process_2023,   "GO BP - Concordant UP",   "results/figures/GO_UP.png")
plot_enrichment(enrich_down$GO_Biological_Process_2023, "GO BP - Concordant DOWN", "results/figures/GO_DOWN.png")
plot_enrichment(enrich_up$KEGG_2021_Human,   "KEGG - Concordant UP",   "results/figures/KEGG_UP.png")
plot_enrichment(enrich_down$KEGG_2021_Human, "KEGG - Concordant DOWN", "results/figures/KEGG_DOWN.png")

cat("GO UP significant terms:",   nrow(enrich_up$GO_Biological_Process_2023[enrich_up$GO_Biological_Process_2023$Adjusted.P.value < 0.05,]), "\n")
cat("GO DOWN significant terms:", nrow(enrich_down$GO_Biological_Process_2023[enrich_down$GO_Biological_Process_2023$Adjusted.P.value < 0.05,]), "\n")
cat("KEGG UP significant terms:",   nrow(enrich_up$KEGG_2021_Human[enrich_up$KEGG_2021_Human$Adjusted.P.value < 0.05,]), "\n")
cat("KEGG DOWN significant terms:", nrow(enrich_down$KEGG_2021_Human[enrich_down$KEGG_2021_Human$Adjusted.P.value < 0.05,]), "\n")
