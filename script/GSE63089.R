if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
rm(list = ls())
BiocManager::install("org.Hs.eg.db")

library(GEOquery)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)
gset <- getGEO("GSE63089", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]
expr <- exprs(gset)        # expression matrix
pdata <- pData(gset)       # sample metadata
head(expr)
dim(expr)
head(pdata)
group <- factor(ifelse(
  grepl("cancer", pdata$source_name_ch1, ignore.case = TRUE),
  "cancer", "normal"
))
table(group)
design <- model.matrix(~ group)
colnames(design)

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

deg_sig$GeneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(deg_sig),
  column = "SYMBOL",
  keytype = "PROBEID",   # change if ENTREZID
  multiVals = "first"
)

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

# Export to CSV
write.csv(deg_final, "DEG_Cancer_vs_Normal_89.csv", row.names = TRUE)
write.csv(upregulated, "Up_DEG_Cancer_vs_Normal_89.csv", row.names = TRUE)
write.csv(downregulated, "Down_DEG_Cancer_vs_Normal_89.csv", row.names = TRUE)









