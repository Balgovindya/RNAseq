#bubbleplot for non model  species using go term enrichment data from David database
library(ggplot2)
# directory
dir_path <- "D:/rnaseq/bal/diffdata/GO/0vs1/"

files <- list.files(dir_path, pattern = "*.txt", full.names = TRUE)

# Read all files and combine into a single data frame
go_data <- do.call(rbind, lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t")
}))
library(dplyr)

go_summary <- go_data %>%
  group_by(Category) %>%
  summarize(Count = n(), .groups = "drop")
library(ggplot2)

ggplot(go_summary, aes(x = reorder(Category, -Count), y = Count)) +
  geom_point(size = 4) +
  coord_flip() +
  labs(title = "GO Term Dot Plot", x = "GO Term", y = "Count") +
  theme_minimal()
top_go_terms <- go_data %>%
  group_by(Category) %>%
  top_n(10, Count) %>%
  ungroup() %>%
  arrange(Category, desc(Count))
ggplot(top_go_terms, aes(x = reorder(Term, -Count), y = Count, color = Category)) +
  geom_point(size = 4) +
  facet_grid(~ Category, scales = "free_y") +
  coord_flip() +
  labs(title = "Top 10 GO Terms by Category", x = "GO Term", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")
#color based on pvalue
ggplot(top_go_terms, aes(x = reorder(Term, -Count), y = Count, color = PValue, size = Count)) +
  geom_point(alpha = 0.7) +
  facet_grid(. ~ Category) +
  coord_flip() +
  labs(title = "1h M treatment", x = "GO Term", y = "Count") +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(2, 6)) +
  theme_minimal() +
  theme(legend.title = element_text("P-Value"),
        axis.text.y = element_text(size = 8))
write.csv(top_go_terms,"D:/rnaseq/bal/diffdata/GO/0vs1/goterms1h.csv")
ggsave("top_GO_terms_1h M treatment.png", width = 12, height = 8)

#0vs3h
dir_path <- "D:/rnaseq/bal/diffdata/GO/0vs3"

files <- list.files(dir_path, pattern = "*.txt", full.names = TRUE)

# Read all files and combine into a single data frame
go_data <- do.call(rbind, lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t")
}))
library(dplyr)

go_summary <- go_data %>%
  group_by(Category) %>%
  summarize(Count = n(), .groups = "drop")
library(ggplot2)


top_go_terms3h <- go_data %>%
  group_by(Category) %>%
  top_n(10, Count) %>%
  ungroup() %>%
  arrange(Category, desc(Count))

#color based on pvalue
ggplot(top_go_terms3h, aes(x = reorder(Term, -Count), y = Count, color = PValue, size = Count)) +
  geom_point(alpha = 0.7) +
  facet_grid(. ~ Category) +
  coord_flip() +
  labs(title = "3h M treatment", x = "GO Term", y = "Count") +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(2, 6)) +
  theme_minimal() +
  theme(legend.title = element_text("P-Value"),
        axis.text.y = element_text(size = 8))
write.csv(top_go_terms3h,"D:/rnaseq/bal/diffdata/GO/0vs3/goterms3h.csv")
ggsave("top_GO_terms_3h M treatment.png", width = 12, height = 8)
#0vs5h
dir_path <- "D:/rnaseq/bal/diffdata/GO/0vs5"

files <- list.files(dir_path, pattern = "*.txt", full.names = TRUE)

# Read all files and combine into a single data frame
go_data <- do.call(rbind, lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t")
}))
library(dplyr)

go_summary <- go_data %>%
  group_by(Category) %>%
  summarize(Count = n(), .groups = "drop")
library(ggplot2)


top_go_terms5h <- go_data %>%
  group_by(Category) %>%
  top_n(10, Count) %>%
  ungroup() %>%
  arrange(Category, desc(Count))
write.csv(top_go_terms5h,"D:/rnaseq/bal/diffdata/GO/0vs5/goterms5h.csv")
#color based on pvalue
ggplot(top_go_terms, aes(x = reorder(Term, -Count), y = Count, color = PValue, size = Count)) +
  geom_point(alpha = 0.7) +
  facet_grid(. ~ Category) +
  coord_flip() +
  labs(title = "5h M treatment", x = "GO Term", y = "Count") +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(2, 6)) +
  theme_minimal() +
  theme(legend.title = element_text("P-Value"),
        axis.text.y = element_text(size = 8))
ggsave("top_GO_terms_5h M treatment.png", width = 12, height = 8)
#########goplot PME data######

setwd ("D:/rnaseq/bal/diffdata/GOterm")
getwd()
# List all text files in the directory
files <- list.files(dir_path, pattern = "*.txt", full.names = TRUE)

# Read all files and combine into a single data frame
go_data <- do.call(rbind, lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t")
}))
library(dplyr)

go_summary <- go_data %>%
  group_by(Category) %>%
  summarize(Count = n(), .groups = "drop")
library(ggplot2)


top_go_terms <- go_data %>%
  group_by(Category) %>%
  top_n(10, Count) %>%
  ungroup() %>%
  arrange(Category, desc(Count))

#color based on pvalue
ggplot(top_go_terms, aes(x = reorder(Term, -Count), y = Count, color = PValue, size = Count)) +
  geom_point(alpha = 0.7) +
  facet_grid(. ~ Category) +
  coord_flip() +
  labs(title = "GOterm PME", x = "GO Term", y = "Count") +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size(range = c(2, 6)) +
  theme_minimal() +
  theme(legend.title = element_text("P-Value"),
        axis.text.y = element_text(size = 8))
ggsave("GOterm PME", width = 12, height = 8)
# Specify the GO term you want to extract genes from
cell_wall_modi <- "GO:0042545~cell wall modification"

# Filter for the specific GO term
specific_go_term <- top_go_terms %>% filter(Term == cell_wall_modi )

# Extract the 'Genes' column for the filtered GO term
genes_list <- specific_go_term %>% pull(Genes) %>% strsplit("/")

# Print the list of genes
print(genes_list)
# Create a volcano plot
results_df=read.csv("D:/rnaseq/bal/diffdata/deseq2_resultsPME.csv")
ggplot(results_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.4, size = 2) +
  theme_minimal() +
  labs(title = "Volcano Plot PME", x = "log2FoldChange", y = "-Log10(padj)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")
# Add a column for significance
results_df$significant <- ifelse(results_df$pvalue < 0.05 & abs(results_df$log2FoldChange) > 1, "Significant", "Non-significant")
library(ggplot2)

# Create volcano plot with color based on significance
# Add a column to classify the genes
results_df$regulation <- ifelse(results_df$padj < 0.05 & results_df$log2FoldChange > 1, "Upregulated",
                                ifelse(results_df$padj < 0.05 & results_df$log2FoldChange < -1, "Downregulated", "Non-significant"))
ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +http://127.0.0.1:35351/graphics/plot_zoom_png?width=1888&height=931
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(P-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")
# Select the top 10 most significant genes by p-value
# Select the top 10 most significant genes by p-value
top_genes <- results_df %>%
  arrange(padj) %>%
  head(10)


# Create the volcano plot with labels for the top 10 genes
library(ggplot2)

# Create the volcano plot with labels for the top 10 genes (GeneID)
ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot PME", x = "Log2 Fold Change", y = "-Log10(padj)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  # Add text labels for top 10 genes by GeneID
  geom_text(data = top_genes, aes(label = GeneID), vjust = 1.5, color = "black", size = 3)
