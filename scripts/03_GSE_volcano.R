library(ggplot2)

res <- read.csv("results/tables/GSE19804_full_DEG_table.csv")

res$Significant <- ifelse(res$adj.P.Val < 0.05 & abs(res$logFC) > 1, "Yes", "No")

p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot - GSE19804",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")

ggsave("results/figures/GSE19804_volcano.png", p, width = 6, height = 5)

cat("Volcano plot saved\n")
