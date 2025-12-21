if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
rm(list = ls())
BiocManager::install("org.Hs.eg.db")

# required library
library(GEOquery)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)

#load data
gset <- getGEO("GSE63089", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]

expr <- exprs(gset)        # expression matrix
pdata <- pData(gset)       # sample metadata

head(expr)
dim(expr)
head(pdata)

#make group cancer vs normal
group <- factor(ifelse(
  grepl("cancer", pdata$source_name_ch1, ignore.case = TRUE),
  "cancer", "normal"
))
table(group)

#make design  matrix
design <- model.matrix(~ group)
colnames(design)

#fit linear model
fit <- lmFit(expr, design)
fit <- eBayes(fit)
deg <- topTable(
  fit,
  coef = "groupnormal",
  number = Inf,
  adjust.method = "BH"
)

head(deg)
deg$logFC_Cancer_vs_Normal <- -deg$logFC

deg_sig <- subset(
  deg,
  adj.P.Val < 0.05 & abs(logFC_Cancer_vs_Normal) > 1
)

upregulated <- subset(deg_sig, logFC_Cancer_vs_Normal > 1)
downregulated <- subset(deg_sig, logFC_Cancer_vs_Normal < -1)

cat("Total DEGs:", nrow(deg_sig), "\n")
cat("Upregulated in cancer:", nrow(upregulated), "\n")
cat("Downregulated in cancer:", nrow(downregulated), "\n")

rownames(expr)[1]
keytypes(org.Hs.eg.db)
fdata <- fData(gset)
head(fdata)
colnames(fdata)

grep("symbol|gene", colnames(fdata), ignore.case = TRUE, value = TRUE)

# Extract only the gene symbol (middle part between the first and second '//')
deg_sig$GeneName <- sapply(
  strsplit(fdata[rownames(deg_sig), "gene_assignment"], " // "),
  function(x) x[2]  # pick the second element
)

head(deg_sig[, c("GeneName", "logFC_Cancer_vs_Normal", "adj.P.Val")])
colnames(deg_sig)

# Keep only the desired columns
deg_final <- deg_sig[, c("GeneName","logFC", "AveExpr"  ,"t","P.Value", "adj.P.Val","logFC_Cancer_vs_Normal",  "B"  )]
upregulated <- subset(deg_final, logFC_Cancer_vs_Normal > 1)
downregulated <- subset(deg_final, logFC_Cancer_vs_Normal < -1)



# Check
head(deg_final)
sum(is.na(deg_final$GeneSymbol))

# Export to CSV
write.csv(deg_final, "F:/DATA_RNASeq/Cancer-Biology/data/GSE63089/DEG_Cancer_vs_Normal_89.csv", row.names = TRUE)
write.csv(upregulated, "F:/DATA_RNASeq/Cancer-Biology/data/GSE63089/Up_DEG_Cancer_vs_Normal_89.csv", row.names = TRUE)
write.csv(downregulated, "F:/DATA_RNASeq/Cancer-Biology/data/GSE63089/Down_DEG_Cancer_vs_Normal_89.csv", row.names = TRUE)

#Top 100 gene

top100_deg <- deg_final[order(deg_final$adj.P.Val), ][1:100, ]

top100_up <- deg_final[
  deg_final$logFC > 1,
][order(deg_final$adj.P.Val[deg_final$logFC > 1]), ][1:100, ]

write.csv(top100_deg,  "F:/DATA_RNASeq/Cancer-Biology/data/GSE63089/Top100_DEGs_89.csv")
write.csv(top100_up,   "F:/DATA_RNASeq/Cancer-Biology/data/GSE63089/Top100_Upregulated_89.csv")










