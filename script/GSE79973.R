if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
rm(list = ls())

BiocManager::install("org.Hs.eg.db")

library(GEOquery)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)

gset <- getGEO("GSE79973", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]
expr <- exprs(gset)        # expression matrix
pdata <- pData(gset)       # sample metadata

head(expr)
dim(expr)
head(pdata)

group <- factor(ifelse(
  grepl("adenocarcinoma", pdata$source_name_ch1, ignore.case = TRUE),
  "adenocarcinoma", "mucosa"
))

group <- factor(group,
                levels = c("mucosa", "adenocarcinoma"),
                labels = c("normal", "tumor"))

table(group)
design <- model.matrix(~ group)
colnames(design)

fit <- lmFit(expr, design)
fit <- eBayes(fit)
deg <- topTable(
  fit,
  coef = "grouptumor",
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

annotation(gset)
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)

head(rownames(deg_sig))

deg_sig$GeneSymbol <- mapIds(
  hgu133plus2.db,
  keys = rownames(deg_sig),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)
deg_sig$GeneName <- mapIds(
  hgu133plus2.db,
  keys = rownames(deg_sig),
  column = "GENENAME",
  keytype = "PROBEID",
  multiVals = "first"
)
deg_sig$ENTREZID <- mapIds(
  hgu133plus2.db,
  keys = rownames(deg_sig),
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

colnames(deg_sig)

# Keep only the desired columns
deg_final <- deg_sig[, c("GeneSymbol","logFC", "AveExpr"  ,"t","P.Value", "adj.P.Val","logFC_Cancer_vs_Normal",  "B" , "GeneName", "ENTREZID")]
sum(is.na(deg_final$GeneSymbol))

deg_final11 <- deg_final[!is.na(deg_final$GeneSymbol), ]
sum(is.na(deg_final11$GeneSymbol))

upregulated <- subset(deg_final11, logFC_Cancer_vs_Normal > 1)
downregulated <- subset(deg_final11, logFC_Cancer_vs_Normal < -1)

# Check
head(deg_final11)

# Export to CSV
write.csv(
  deg_final11,
  file = "F:/DATA_RNASeq/Cancer-Biology/data/GSE79973/DEG_Cancer_vs_Normal_73.csv",
  row.names = TRUE
)
write.csv(
  upregulated, 
  file = "F:/DATA_RNASeq/Cancer-Biology/data/GSE79973/Up_DEG_Cancer_vs_Normal_73.csv", 
  row.names = TRUE
)
write.csv(
  downregulated, 
  file= "F:/DATA_RNASeq/Cancer-Biology/data/GSE79973/Down_DEG_Cancer_vs_Normal_73.csv", 
  row.names = TRUE
)

# Top 100 gene
top100_deg <- deg_final11[order(deg_final11$adj.P.Val), ][1:100, ]

top100_up <- deg_final11[
  deg_final11$logFC > 1,
][order(deg_final11$adj.P.Val[deg_final11$logFC > 1]), ][1:100, ]

write.csv(top100_deg,  "F:/DATA_RNASeq/Cancer-Biology/data/GSE79973/Top100_DEGs_73.csv")
write.csv(top100_up,   "F:/DATA_RNASeq/Cancer-Biology/data/GSE79973/Top100_Upregulated_73.csv")








