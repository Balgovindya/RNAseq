library(pheatmap)
library(reshape2)

# Define the list of target GO terms
target_go_terms <- c("GO:0006355", "GO:0008017", "GO:0051301", "GO:0006468", "GO:0005975")

file_1h <- read.csv('D:/rnaseq/bal/diffdata/GO/0vs1/goterm_zscore_1h.csv', header = TRUE, row.names = 1)
file_3h <- read.csv('D:/rnaseq/bal/diffdata/GO/0vs3/goterm_zscore_3h.csv', header = TRUE, row.names = 1)
file_5h <- read.csv('D:/rnaseq/bal/diffdata/GO/0vs5/goterm_zscore_5h.csv', header = TRUE, row.names = 1)

# Clean data
file_1h <- na.omit(file_1h)
file_3h <- na.omit(file_3h)
file_5h <- na.omit(file_5h)

# indicate time points
colnames(file_1h)[colnames(file_1h) == "logFC"] <- "log2FoldChange_1h"
colnames(file_3h)[colnames(file_3h) == "logFC"] <- "log2FoldChange_3h"
colnames(file_5h)[colnames(file_5h) == "logFC"] <- "log2FoldChange_5h"

# Merge data on gene names
combined_data <- merge(file_1h, file_3h, by = "genes", all = TRUE)
combined_data <- merge(combined_data, file_5h, by = "genes", all = TRUE)
combined_data <- na.omit(combined_data)

# Loop through each GO term for heatmap
for (go_term in target_go_terms) {
  go_data <- subset(combined_data, ID == go_term, 
                    select = c(genes, log2FoldChange_1h, log2FoldChange_3h, log2FoldChange_5h))
  
  if (nrow(go_data) > 0) {
    # Prepare data for heatmap (exclude 'genes' column)
    heatmap_data <- as.matrix(go_data[, -1])
    rownames(heatmap_data) <- go_data$genes
    
    # Save each heatmap as a PDF
    pdf_filename <- paste0("Heatmap_", go_term, ".pdf")
    pdf(pdf_filename)
    
    # Generate the heatmap
    pheatmap(heatmap_data,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             scale = "none",
             color = colorRampPalette(c("blue", "white", "purple"))(50),
             main = paste("Heatmap of", go_term, "log2FoldChange at Different Time Points"),
             fontsize_row = 10,
             fontsize_col = 12,
             cellheight = 20)
    
    dev.off()
  } else {
    message(paste("No data available for", go_term))
  }
}
