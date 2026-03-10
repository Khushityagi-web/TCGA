# 19_PPI_network.R
library(STRINGdb)
library(ggplot2)

res   <- read.csv("results/tables/Survival_screen_all.csv")
sig85 <- res[res$FDR < 0.05, ]
genes <- sig85$Gene
cat("Building PPI for", length(genes), "genes\n")

string_db <- STRINGdb$new(version="11.5", species=9606,
                           score_threshold=400, input_directory="data_raw/")

mapped <- string_db$map(data.frame(gene=genes), "gene", removeUnmappedRows=TRUE)
cat("Mapped to STRING:", nrow(mapped), "genes\n")

interactions <- string_db$get_interactions(mapped$STRING_id)
write.csv(interactions, "results/tables/PPI_interactions.csv", row.names=FALSE)

# Network plot
png("results/figures/PPI_network.png", width=1200, height=1000)
string_db$plot_network(mapped$STRING_id)
dev.off()

# Hub genes by degree
all_nodes <- c(interactions$from, interactions$to)
degree    <- sort(table(all_nodes), decreasing=TRUE)
hub_df    <- data.frame(STRING_id=names(degree), Degree=as.numeric(degree))
hub_df    <- merge(hub_df, mapped[, c("gene","STRING_id")], by="STRING_id")
hub_df    <- hub_df[order(-hub_df$Degree), ]

write.csv(hub_df, "results/tables/Hub_genes.csv", row.names=FALSE)
cat("\nTop 10 hub genes:\n")
print(head(hub_df[, c("gene","Degree")], 10))
