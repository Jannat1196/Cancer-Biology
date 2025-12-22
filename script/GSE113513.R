if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
rm(list = ls())
BiocManager::install("hugene10sttranscriptcluster.db")

library(hugene10sttranscriptcluster.db)
library(GEOquery)
library(limma)
library(AnnotationDbi)

#Load the Data 
gset <- getGEO("GSE113513", GSEMatrix = TRUE)
gset <- gset[[1]]

#For Gene name adding need platform
annotation(gset)
# Should return "GPL15207"
gpl <- getGEO("GPL15207", destdir=".")
gpl_table <- Table(gpl)
head(gpl_table)

expr <- exprs(gset)        # expression matrix
pdata <- pData(gset)       # sample metadata

head(expr)
dim(expr)
head(pdata)
colnames(pdata)

# Group data Tumor vs Normal
group <- factor(ifelse(
  pdata$`tissue:ch1` == "colorectal carcinoma",
  "colorectal carcinoma",
  "non-cancerous colorectal"
))

group <- factor(group,
                levels = c("non-cancerous colorectal", "colorectal carcinoma"),
                labels = c("normal", "tumor"))

table(group)

#create matrix
design <- model.matrix(~ group)
colnames(design)

# fit the data
fit <- lmFit(expr, design)
fit <- eBayes(fit)

#Find DEG
deg <- topTable(
  fit,
  coef = "grouptumor",
  number = Inf,
  adjust.method = "BH"
)

head(deg)
deg$logFC_Tumor_vs_Normal <- -deg$logFC
# Find Significant DEG at 5% level of sig. and logcFold >1 
deg_sig <- subset(
  deg,
  adj.P.Val < 0.05 & abs(logFC_Tumor_vs_Normal) > 1
)

upregulated <- subset(deg_sig, logFC_Tumor_vs_Normal > 1)
downregulated <- subset(deg_sig, logFC_Tumor_vs_Normal < -1)

#summary 
cat("Total DEGs:", nrow(deg_sig), "\n")
cat("Upregulated in cancer:", nrow(upregulated), "\n")
cat("Downregulated in cancer:", nrow(downregulated), "\n")

# Add gene name to merge with deg_sig

# Suppose your DEGs table is deg_sig with rownames = probe IDs
deg_sig$ProbeID <- rownames(deg_sig)

# Merge with platform table
deg_sig_annotated <- merge(
  deg_sig,
  gpl_table[, c("ID", "Gene Symbol", "Gene Title", "Entrez Gene")],
  by.x = "ProbeID",
  by.y = "ID",
  all.x = TRUE
)
colnames(gpl_table)

head(deg_sig_annotated[, c("Gene Symbol", "Gene Title" , "Entrez Gene")])

# Rename columns if needed
colnames(deg_sig_annotated)[(ncol(deg_sig)+1):ncol(deg_sig_annotated)] <- c("Gene_Symbol", "Gene_Name", "ENTREZ_ID")

# Keep only the desired columns
deg_final <- deg_sig_annotated[, c( "ProbeID", "Gene_Symbol", "logFC","AveExpr" ,"t", "P.Value","adj.P.Val", "B", "logFC_Tumor_vs_Normal",       
                          "Gene_Name" , "ENTREZ_ID"  )]

#Find any NA value
sum(is.na(deg_final$GeneSymbol))

upregulated <- subset(deg_final, logFC_Tumor_vs_Normal > 1)
downregulated <- subset(deg_final, logFC_Tumor_vs_Normal < -1)

#summary
nrow(deg_final)
nrow( upregulated)
nrow(downregulated)

# Check
head(deg_final)

# Export to CSV
write.csv(
  deg_final,
  file = "F:/DATA_RNASeq/Cancer-Biology/data/GSE113513/DEG_Cancer_vs_Normal_13.csv",
  row.names = TRUE
)
write.csv(
  upregulated, 
  file = "F:/DATA_RNASeq/Cancer-Biology/data/GSE113513/Up_DEG_Cancer_vs_Normal_13.csv", 
  row.names = TRUE
)
write.csv(
  downregulated, 
  file= "F:/DATA_RNASeq/Cancer-Biology/data/GSE113513/Down_DEG_Cancer_vs_Normal_13.csv", 
  row.names = TRUE
)

# Find 100 gene
top100_deg <- deg_final[order(deg_final$adj.P.Val), ][1:100, ]

#Find Top 100 up regulated gene
top100_up <- deg_final[
  deg_final$logFC > 1,
][order(deg_final$adj.P.Val[deg_final$logFC > 1]), ][1:100, ]

#Export  the top gene
write.csv(top100_deg,  "F:/DATA_RNASeq/Cancer-Biology/data/GSE113513/Top100_DEGs_13.csv")
write.csv(top100_up,   "F:/DATA_RNASeq/Cancer-Biology/data/GSE113513/Top100_Upregulated_13.csv")


#file.info("script/GSE44861.R")$size

