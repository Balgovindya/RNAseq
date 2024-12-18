#!/usr/bin/env Rscript
setwd("D:/rnaseq/bal")
# Load necessary libraries
library("DESeq2")


# input files
countData <-read.csv("data.csv", header = T)
samples_to_extract <- c("Geneid", "WTM01", "WTM02", "WTM03", "WTM51", "WTM52", "WTM53")

# keep only the columns of interest
countData_subset5 <- countData[, colnames(countData) %in% samples_to_extract]

countData5<-data.matrix(countData_subset5[,2:ncol(countData_subset5)])
rnames<-countData_subset5[,1]
rownames(countData5) <-rnames
# Convert the count data (excluding the first column) into a matrix
countData5 <- data.matrix(countData_subset5[, 2:ncol(countData_subset5)])

# Extract the gene IDs
rnames <- countData_subset5[, 1]
rownames(countData5) <- rnames

# Now countData5 is a matrix with gene IDs as row names

# Assign a name to the first column
colnames(df_set1)[1] <- "Gene_ID"

# to check if the column name has been assigned
colnames(df_set1)


colnames(countData5)[0] <- "GENE_ID"
colData5 <- data.frame(
  id = c("WTM01", "WTM02", "WTM03", "WTM51", "WTM52", "WTM53"),
  dex = c("control", "control", "control", "treatment3", "treatment3", "treatment3"))
rownames(colData5) <- colData5$id
countData<-data.matrix(countData[,2:ncol(countData)])
rnames<-countData[,1]
rownames(countData) <-rnames

colData <- data.frame(
  id = c("WTM01", "WTM02", "WTM03", "WTM11", "WTM12", "WTM13", 
         "WTM31", "WTM32", "WTM33", "WTM51", "WTM52", "WTM53"),
  dex = c("control", "control", "control", 
          "treatment1", "treatment1", "treatment1", 
          "treatment2", "treatment2", "treatment2", 
          "treatment3", "treatment3", "treatment3")
)
rownames(colData) <- colData$id
# Step 1: Load the data
all(rownames(colData) == colnames(countData))

# Ensure sample names in the metadata match the column names in the count matrix
all(rownames(colData5) == colnames(countData5))  

# Step 2: Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ dex)

# Step 3: Filter low-count genes (optional but recommended)
# Keep genes with at least 10 reads across all samples
dds <- dds[rowSums(counts(dds)) > 10, ]

# Step 4: Run the DESeq2 pipeline
dds5 <- DESeq(dds5)
dds <- DESeq(dds)
# Step 5: Results extraction
# Perform the differential expression analysis between "treatment" and "control"
res5 <- results(dds5)
res <- results(dds)

# Step 6: View and export results
# Sort results by adjusted p-value
resOrdered5 <- res5[order(res5$padj), ]
resOrdered <- res[order(res$padj), ]
# View summary of results
summary(res5)
summary(res)
# Filter for upregulated genes (log2FoldChange > 1)
upregulated_genes5 <- res5[res5$log2FoldChange > 2 & res5$padj < 0.05, ]

# Filter for downregulated genes (log2FoldChange < -1)
downregulated_genes5 <- res5[res5$log2FoldChange < -2 & res5$padj < 0.05, ]
diff5=merge(upregulated_genes5,downregulated_genes5)
# Write results to a CSV file
write.csv(as.data.frame(resOrdered5), file = "deseq2_results5h.csv")

# Step 7: Plot MA-plot (log fold change vs mean expression)
plotMA(res5, main="5h M treatment", ylim=c(-2, 2))
plotMA(res, main="M treatment", ylim=c(-2, 2))

# Step 8: Plot normalized counts for a specific gene (optional)
# You can plot the normalized counts for a gene of interest
# Example: Plotting normalized counts for "GENE_ID"
plotCounts(dds5, gene="", intgroup="condition")
df_set1
plotCounts(dds, gene="", intgroup="condition")
# Step 9: PCA plot for sample clustering
# Transform the data using regularized-log transformation (rlog)
rld <- rlog(dds5, blind = TRUE)
rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup="id")
ggsave("PCAall.png", width = 8, height = 8)



#differential analysis between 0 and 1 hour

# Paths to input files
countData <-read.csv("data.csv", header = T)
samples_to_extract1 <- c("Geneid", "WTM01", "WTM02", "WTM03", "WTM11", "WTM12", "WTM13")

# Subset the countData matrix to keep only the columns of interest
countData_subset1 <- countData[, colnames(countData) %in% samples_to_extract1]

countData1<-data.matrix(countData_subset1[,2:ncol(countData_subset1)])
rnames<-countData_subset1[,1]
rownames(countData1) <-rnames
colData1 <- data.frame(
  id = c("WTM01", "WTM02", "WTM03", "WTM11", "WTM12", "WTM13"),
  dex = c("control", "control", "control", "treatment1", "treatment1", "treatment1"))
rownames(colData1) <- colData1$id


colData <- data.frame(
  id = c("WTM01", "WTM02", "WTM03", "WTM11", "WTM12", "WTM13", 
         "WTM31", "WTM32", "WTM33", "WTM51", "WTM52", "WTM53"),
  dex = c("control", "control", "control", 
          "treatment1", "treatment1", "treatment1", 
          "treatment2", "treatment2", "treatment2", 
          "treatment3", "treatment3", "treatment3")
)
rownames(colData) <- colData$id
# Step 1: Load the data


# Ensure sample names in the metadata match the column names in the count matrix
all(rownames(colData1) == colnames(countData1))  

# Step 2: Create DESeq2 dataset
dds1 <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ dex)

# Step 3: Filter low-count genes (optional but recommended)
# Keep genes with at least 10 reads across all samples
dds1 <- dds1[rowSums(counts(dds1)) > 10, ]

# Step 4: Run the DESeq2 pipeline
dds1 <- DESeq(dds1)

# Step 5: Results extraction
# Perform the differential expression analysis between "treatment" and "control"
res1 <- results(dds1)

# Step 6: View and export results
# Sort results by adjusted p-value
resOrdered1 <- res1[order(res1$padj), ]

# View summary of results
summary(res1)

# Write results to a CSV file
write.csv(as.data.frame(resOrdered1), file = "deseq2_results1h.csv")

# Step 7: Plot MA-plot (log fold change vs mean expression)
plotMA(res1, main="1h M treatment", ylim=c(-2, 2))

# Step 8: Plot normalized counts for a specific gene (optional)
# You can plot the normalized counts for a gene of interest
# Example: Plotting normalized counts for "GENE_ID"
plotCounts(dds1, gene="", intgroup="condition")

# Step 9: PCA plot for sample clustering
# Transform the data using regularized-log transformation (rlog)
rld <- rlog(dds1, blind = TRUE)

# PCA plot of samples
plotPCA(rld, intgroup="id")
ggsave("PCA1h.png", width = 8, height = 8)
#differential analysis between 0 and 3 hour



# Paths to input files
countData <-read.csv("data.csv", header = T)  
# Define the samples you want to extract
samples_to_extract3 <- c("Geneid", "WTM01", "WTM02", "WTM03", "WTM31", "WTM32", "WTM33")

# Subset the countData matrix to keep only the columns of interest
countData_subset3 <- countData[, colnames(countData) %in% samples_to_extract3]

countData3<-data.matrix(countData_subset3[,2:ncol(countData_subset3)])
rnames<-countData_subset3[,1]
rownames(countData3) <-rnames
colData3 <- data.frame(
  id = c("WTM01", "WTM02", "WTM03", "WTM31", "WTM32", "WTM33"),
  dex = c("control", "control", "control", "treatment1", "treatment2", "treatment2"))
rownames(colData3) <- colData3$id


colData <- data.frame(
  id = c("WTM01", "WTM02", "WTM03", "WTM11", "WTM12", "WTM13", 
         "WTM31", "WTM32", "WTM33", "WTM51", "WTM52", "WTM53"),
  dex = c("control", "control", "control", 
          "treatment1", "treatment1", "treatment1", 
          "treatment2", "treatment2", "treatment2", 
          "treatment3", "treatment3", "treatment3")
)
rownames(colData) <- colData$id
# Step 1: Load the data


# Ensure sample names in the metadata match the column names in the count matrix
all(rownames(colData3) == colnames(countData3))  

# Step 2: Create DESeq2 dataset
dds3 <- DESeqDataSetFromMatrix(countData = countData3, colData = colData3, design = ~ dex)

# Step 3: Filter low-count genes (optional but recommended)
# Keep genes with at least 10 reads across all samples
dds3 <- dds3[rowSums(counts(dds3)) > 10, ]

# Step 4: Run the DESeq2 pipeline
dds3 <- DESeq(dds3)

# Step 5: Results extraction
# Perform the differential expression analysis between "treatment3" and "control"
res3 <- results(dds3)

# Step 6: View and export results
# Sort results by adjusted p-value
resOrdered3 <- res3[order(res3$padj), ]

# View summary of results
summary(res3)

# Write results to a CSV file
write.csv(as.data.frame(resOrdered3), file = "deseq2_results3h.csv")

# Step 7: Plot MA-plot (log fold change vs mean expression)
plotMA(res3, main="3h M treatment", ylim=c(-2, 2))

# Step 8: Plot normalized counts for a specific gene (optional)
# You can plot the normalized counts for a gene of interest
# Example: Plotting normalized counts for "GENE_ID"
plotCounts(dds3, gene="", intgroup="condition")

# Step 9: PCA plot for sample clustering
# Transform the data using regularized-log transformation (rlog)
rld3 <- rlog(dds3, blind = TRUE)

# PCA plot of samples
plotPCA(rld3, intgroup="id")

ggsave("PCA3h.png", width = 8, height = 8)

# Paths to input files
#differential analysis between 1 and 3 hour

countData <-read.csv("data.csv", header = T)  
# Define the samples you want to extract
samples_to_extract13 <- c("Geneid", "WTM11", "WTM12", "WTM13", "WTM31", "WTM32", "WTM33")

# Subset the countData matrix to keep only the columns of interest
countData_subset13 <- countData[, colnames(countData) %in% samples_to_extract13]

countData13<-data.matrix(countData_subset13[,2:ncol(countData_subset13)])
rnames<-countData_subset13[,1]
rownames(countData13) <-rnames
colData13 <- data.frame(
  id = c("WTM11", "WTM12", "WTM13", "WTM31", "WTM32", "WTM33"),
  dex = c("control", "control", "control", "treatment1", "treatment2", "treatment2"))
rownames(colData13) <- colData13$id




# Ensure sample names in the metadata match the column names in the count matrix
all(rownames(colData13) == colnames(countData13))  

# Step 2: Create DESeq2 dataset
dds13 <- DESeqDataSetFromMatrix(countData = countData13, colData = colData13, design = ~ dex)

# Step 3: Filter low-count genes (optional but recommended)
# Keep genes with at least 10 reads across all samples
dds13 <- dds13[rowSums(counts(dds13)) > 10, ]

# Step 4: Run the DESeq2 pipeline
dds13 <- DESeq(dds13)

# Step 5: Results extraction
# Perform the differential expression analysis between "treatment" and "control"
res13 <- results(dds13)

# Step 6: View and export results
# Sort results by adjusted p-value
resOrdered13 <- res13[order(res13$padj), ]

# View summary of results
summary(res13)

# Write results to a CSV file
write.csv(as.data.frame(resOrdered13), file = "deseq2_results13h.csv")

# Step 7: Plot MA-plot (log fold change vs mean expression)
plotMA(res3, main="1&3h M treatment", ylim=c(-2, 2))

# Step 8: Plot normalized counts for a specific gene (optional)
# You can plot the normalized counts for a gene of interest
# Example: Plotting normalized counts for "GENE_ID"
plotCounts(dds3, gene="", intgroup="condition")

# Step 9: PCA plot for sample clustering
# Transform the data using regularized-log transformation (rlog)
rld13 <- rlog(dds13, blind = TRUE)

# PCA plot of samples
plotPCA(rld13, intgroup="id")
ggsave("PCA13h.png", width = 8, height = 8)
#Perform the differential expression analysis between 5 and 3

countData <-read.csv("data.csv", header = T)  
# Define the samples you want to extract
samples_to_extract35 <- c("Geneid", "WTM31", "WTM32", "WTM33", "WTM51", "WTM52", "WTM53")

# Subset the countData matrix to keep only the columns of interest
countData_subset35 <- countData[, colnames(countData) %in% samples_to_extract35]

countData35<-data.matrix(countData_subset35[,2:ncol(countData_subset35)])
rnames<-countData_subset35[,1]
rownames(countData35) <-rnames
colData35 <- data.frame(
  id = c("WTM31", "WTM32", "WTM33", "WTM51", "WTM52", "WTM53"),
  dex = c("control", "control", "control", "treatment2", "treatment2", "treatment2"))
rownames(colData35) <- colData35$id




# Ensure sample names in the metadata match the column names in the count matrix
all(rownames(colData35) == colnames(countData35))  

# Step 2: Create DESeq2 dataset
dds35 <- DESeqDataSetFromMatrix(countData = countData35, colData = colData35, design = ~ dex)

# Step 3: Filter low-count genes (optional but recommended)
# Keep genes with at least 10 reads across all samples
dds35 <- dds35[rowSums(counts(dds35)) > 10, ]

# Step 4: Run the DESeq2 pipeline
dds35 <- DESeq(dds35)

# Step 5: Results extraction
# Perform the differential expression analysis between "treatment" and "control"
res35 <- results(dds35)

# Step 6: View and export results
# Sort results by adjusted p-value
resOrdered35 <- res35[order(res35$padj), ]

# View summary of results
summary(res35)


# Write results to a CSV file
write.csv(as.data.frame(resOrdered35), file = "deseq2_results35h.csv")

# Step 7: Plot MA-plot (log fold change vs mean expression)
plotMA(res35, main="3&5h M treatment", ylim=c(-2, 2))

# Step 8: Plot normalized counts for a specific gene (optional)
# You can plot the normalized counts for a gene of interest
# Example: Plotting normalized counts for "GENE_ID"
plotCounts(dds35, gene="", intgroup="condition")

# Step 9: PCA plot for sample clustering
# Transform the data using regularized-log transformation (rlog)
rld35 <- rlog(dds35, blind = TRUE)

# PCA plot of samples
plotPCA(rld35, intgroup="id")
ggsave("PCA35h.png", width = 8, height = 8)
 #Perform the differential expression analysis between 5 and 1

countData <-read.csv("data.csv", header = T)  
# Define the samples you want to extract
samples_to_extract15 <- c("Geneid", "WTM11", "WTM12", "WTM13", "WTM51", "WTM52", "WTM53")

# Subset the countData matrix to keep only the columns of interest
countData_subset15 <- countData[, colnames(countData) %in% samples_to_extract15]

countData15<-data.matrix(countData_subset15[,2:ncol(countData_subset15)])
rnames<-countData_subset15[,1]
rownames(countData15) <-rnames
colData15 <- data.frame(
  id = c("WTM11", "WTM12", "WTM13", "WTM51", "WTM52", "WTM53"),
  dex = c("control", "control", "control", "treatment2", "treatment2", "treatment2"))
rownames(colData15) <- colData15$id




# Ensure sample names in the metadata match the column names in the count matrix
all(rownames(colData15) == colnames(countData15))  

# Step 2: Create DESeq2 dataset
dds15 <- DESeqDataSetFromMatrix(countData = countData15, colData = colData15, design = ~ dex)

# Step 3: Filter low-count genes (optional but recommended)
# Keep genes with at least 10 reads across all samples
dds15 <- dds15[rowSums(counts(dds15)) > 10, ]

# Step 4: Run the DESeq2 pipeline
dds15 <- DESeq(dds15)

# Step 5: Results extraction
# Perform the differential expression analysis between "treatment" and "control"
res15 <- results(dds15)

# Step 6: View and export results
# Sort results by adjusted p-value
resOrdered15 <- res15[order(res15$padj), ]

# View summary of results
summary(res15)

# Write results to a CSV file
write.csv(as.data.frame(resOrdered15), file = "deseq2_results15h.csv")

# Step 7: Plot MA-plot (log fold change vs mean expression)
plotMA(res15, main="1&5h M treatment", ylim=c(-2, 2))

# Step 8: Plot normalized counts for a specific gene (optional)
# You can plot the normalized counts for a gene of interest
# Example: Plotting normalized counts for "GENE_ID"
plotCounts(dds15, gene="", intgroup="condition")

# Step 9: PCA plot for sample clustering
# Transform the data using regularized-log transformation (rlog)
rld15 <- rlog(dds15, blind = TRUE)

# PCA plot of samples
plotPCA(rld15, intgroup="id")
ggsave("PCA35h.png", width = 8, height = 8)
#########################################################################DIFFERENTIAL ANALYSIS  DATA############
# Paths to input files
library("DESeq2")
getwd()
setwd("D:/rnaseq//diffdata")
countData <-read.csv("data.csv", header = T)  
# Define the samples you want to extract
samples_to_extract3 <- c("Geneid", "EV1","EV2","EV3","EV4","PME13","PME1","PME20","PME7"
)

# Subset the countData matrix to keep only the columns of interest
countData_subset3 <- countData[, colnames(countData) %in% samples_to_extract3]

countData3<-data.matrix(countData_subset3[,2:ncol(countData_subset3)])
rnames<-countData_subset3[,1]
rownames(countData3) <-rnames
colData <- data.frame(
  id = c("EV1","EV2","EV3","EV4","PME13","PME1","PME20","PME7"),
  dex = c("control", "control", "control","control", "treatment1","treatment1", "treatment1", "treatment1"))
rownames(colData) <- colData$id




# Step 1: Load the data


# Ensure sample names in the metadata match the column names in the count matrix
all(rownames(colData) == colnames(countData3))  

# Step 2: Create DESeq2 dataset
ddsPME <- DESeqDataSetFromMatrix(countData = countData3, colData = colData, design = ~ dex)

# Step 3: Filter low-count genes (optional but recommended)
# Keep genes with at least 10 reads across all samples
ddsPME <- ddsPME[rowSums(counts(ddsPME)) > 10, ]

# Step 4: Run the DESeq2 pipeline
ddsPME <- DESeq(ddsPME)

# Step 5: Results extraction
# Perform the differential expression analysis between "treatment3" and "control"
resPME <- results(ddsPME)

# Step 6: View and export results
# Sort results by adjusted p-value
resOrderedPME <- resPME[order(resPME$padj), ]

# View summary of results
summary(resPME)

# Write results to a CSV file
write.csv(as.data.frame(resOrderedPME), file = "deseq2_resultsPME.csv")

# Step 7: Plot MA-plot (log fold change vs mean expression)
plotMA(resPME, main="PME", ylim=c(-2, 2))

# Step 8: Plot normalized counts for a specific gene (optional)
# You can plot the normalized counts for a gene of interest
# Example: Plotting normalized counts for "GENE_ID"
plotCounts(ddsPME, gene="", intgroup="condition")

# Step 9: PCA plot for sample clustering
# Transform the data using regularized-log transformation (rlog)
rldPME <- rlog(ddsPME, blind = TRUE)

# PCA plot of samples
plotPCA(rldPME, intgroup="id",main="PME overexpressor PCA")

ggsave("PCA.png", width = 8, height = 8)
filtered_genesPME <- resOrderedPME[resOrderedPME$padj < 0.05 & resOrderedPME$log2FoldChange > -1.5 & resOrderedPME$log2FoldChange < 1.5, ]
dim(filtered_genesPME)
filtered_genesPME <- resOrderedPME[!is.na(resOrderedPME$padj) & resOrderedPME$padj < 0.05 & resOrderedPME$log2FoldChange > -1.5 & resOrderedPME$log2FoldChange < 1.5, ]
