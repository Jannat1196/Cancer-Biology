# ===================== 1 CELL : TCGA COAD + READ RNA-seq =====================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

BiocManager::install("WGCNA")
library(TCGAbiolinks)
library(SummarizedExperiment)  # for assay()
library(edgeR)
library(ggplot2)
library(EnhancedVolcano)  # optional, for volcano plots

get_tcga <- function(project){
  q <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(q)
  GDCprepare(q)
}

# Download and prepare data for both projects
coad <- get_tcga("TCGA-COAD")
read <- get_tcga("TCGA-READ")


# Expression matrices
coad_counts <- assay(coad)
read_counts <- assay(read)

# Summary table
# Use `sample_type` which is more consistent across TCGA projects
coad_type <- colData(coad)$shortLetterCode
if(is.null(coad_type)) coad_type <- colData(coad)$sample_type

read_type <- colData(read)$shortLetterCode
if(is.null(read_type)) read_type <- colData(read)$sample_type

data.frame(
  Dataset = c("TCGA-COAD", "TCGA-READ"),
  No_of_Normal = c(sum(coad_type %in% c("NT", "Solid Tissue Normal")),
                   sum(read_type %in% c("NT", "Solid Tissue Normal"))),
  No_of_Tumor  = c(sum(coad_type %in% c("TP", "Primary Tumor")),
                   sum(read_type %in% c("TP", "Primary Tumor"))),
  PlatformID  = c(if("platform" %in% colnames(colData(coad))) unique(colData(coad)$platform)[1] else NA,
                  if("platform" %in% colnames(colData(read))) unique(colData(read)$platform)[1] else NA),
  No_of_Row   = c(nrow(coad_counts),
                  nrow(read_counts))
)
# Keep genes with counts per million > 1 in at least, say, 10% of samples
keep_coad <- rowSums(cpm(coad_counts) > 1) >= 0.1*ncol(coad_counts)
coad_counts <- coad_counts[keep_coad, ]

keep_read <- rowSums(cpm(read_counts) > 1) >= 0.1*ncol(read_counts)
read_counts <- read_counts[keep_read, ]

# Function to perform edgeR DE analysis
dge_analysis <- function(count_matrix, group_vector){
  group <- factor(
    ifelse(group_vector %in% c("NT","Solid Tissue Normal"), "Normal", "Tumor"), 
    levels = c("Normal","Tumor")
   )
  y <- DGEList(counts = count_matrix, group = group)
  y <- calcNormFactors(y)
  design <- model.matrix(~group)  # Design matrix
  y <- estimateDisp(y, design)  # Estimate dispersion
  fit <- glmQLFit(y, design)   # Fit model
  qlf <- glmQLFTest(fit, coef=2)
  degs <- topTags(qlf, n=Inf)$table
  degs$Significant <- abs(degs$logFC) > 2 & degs$PValue < 0.01
  return(degs)
}

coad_degs <- dge_analysis(coad_counts, coad_type)
read_degs <- dge_analysis(read_counts, read_type)

sum(coad_degs$Significant)
sum(read_degs$Significant)

# Separate Up- and Down-regulated genes
#COAD
coad_up   <- coad_degs[coad_degs$logFC >  2 & coad_degs$PValue < 0.01, ]
coad_down <- coad_degs[coad_degs$logFC < -2 & coad_degs$PValue < 0.01, ]

# READ
read_up   <- read_degs[read_degs$logFC >  2 & read_degs$PValue < 0.01, ]
read_down <- read_degs[read_degs$logFC < -2 & read_degs$PValue < 0.01, ]



### EnhancedVolcano

EnhancedVolcano(coad_degs,
                lab = rownames(coad_degs),
                x = "logFC",
                y = "PValue",
                FCcutoff = 2,
                pCutoff = 0.01,
                title = "TCGA-COAD Volcano Plot")

EnhancedVolcano(read_degs,
                lab = rownames(read_degs),
                x = "logFC",
                y = "PValue",
                FCcutoff = 2,
                pCutoff = 0.01,
                title = "TCGA-READ Volcano Plot")

# COAD
write.csv(coad_degs, "TCGA_COAD_all_DEGs.csv")
write.csv(coad_up,   "TCGA_COAD_upregulated.csv")
write.csv(coad_down, "TCGA_COAD_downregulated.csv")

# READ
write.csv(read_degs, "TCGA_READ_all_DEGs.csv")
write.csv(read_up,   "TCGA_READ_upregulated.csv")
write.csv(read_down, "TCGA_READ_downregulated.csv")

## ===============================
## Add Gene names + Export DEGs
## ===============================

library(biomaRt)

# 1. Connect to Ensembl
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# 2. Function to annotate DEG table
annotate_deg <- function(deg_table){
  
  deg_table$Ensembl_ID <- sub("\\..*", "", rownames(deg_table))
  
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
    filters    = "ensembl_gene_id",
    values     = deg_table$Ensembl_ID,
    mart       = mart
  )
  
  deg_annot <- merge(
    deg_table,
    gene_map,
    by.x = "Ensembl_ID",
    by.y = "ensembl_gene_id",
    all.x = TRUE
  )
  
  return(deg_annot)
}

# 3. Annotate all DEG sets
coad_degs_annot <- annotate_deg(coad_degs)
coad_up_annot   <- annotate_deg(coad_up)
coad_down_annot <- annotate_deg(coad_down)

read_degs_annot <- annotate_deg(read_degs)
read_up_annot   <- annotate_deg(read_up)
read_down_annot <- annotate_deg(read_down)

# 4. Export CSV files
write.csv(coad_degs_annot, "TCGA_COAD_all_DEGs1.csv", row.names = FALSE)
write.csv(coad_up_annot,   "TCGA_COAD_upregulated2.csv", row.names = FALSE)
write.csv(coad_down_annot, "TCGA_COAD_downregulated3.csv", row.names = FALSE)

write.csv(read_degs_annot, "TCGA_READ_all_DEGs1.csv", row.names = FALSE)
write.csv(read_up_annot,   "TCGA_READ_upregulated2.csv", row.names = FALSE)
write.csv(read_down_annot, "TCGA_READ_downregulated3.csv", row.names = FALSE)

# -------------------------------
# Find common DEGs by GeneSymbol
# -------------------------------

# Only significant DEGs (Significant == TRUE)
coad_sig <- coad_degs_annot[coad_degs_annot$Significant == TRUE, ]
read_sig <- read_degs_annot[read_degs_annot$Significant == TRUE, ]

# Get common gene symbols
common_genes <- intersect(
  na.omit(coad_sig$GeneSymbol),
  na.omit(read_sig$GeneSymbol)
)

# Subset DEG tables for common genes
coad_common <- coad_sig[coad_sig$GeneSymbol %in% common_genes, ]
read_common <- read_sig[read_sig$GeneSymbol %in% common_genes, ]


colnames(coad_common)
colnames(read_common)
# Step 1: Select only columns + remove duplicates
coad_common <- coad_sig[coad_sig$hgnc_symbol %in% common_genes, ]
coad_common <- coad_common[, c("hgnc_symbol","logFC","PValue","FDR")]
coad_common <- coad_common[!duplicated(coad_common$hgnc_symbol), ]

read_common <- read_sig[read_sig$hgnc_symbol %in% common_genes, ]
read_common <- read_common[, c("hgnc_symbol","logFC","PValue","FDR")]
read_common <- read_common[!duplicated(read_common$hgnc_symbol), ]

# Step 2: Merge by hgnc_symbol
common_DEGs <- merge(
  coad_common,
  read_common,
  by = "hgnc_symbol",
  suffixes = c("_COAD","_READ")
)

# Step 3: Rename first column
colnames(common_DEGs)[1] <- "GeneSymbol"

# Step 4: Export CSV
write.csv(common_DEGs, "TCGA_COAD_READ_common_DEGs.csv", row.names = FALSE)

# Step 5: Quick check
head(common_DEGs)


# ===================================================
# OFFLINE: Find common DEGs using existing Gene Symbols
# ===================================================

# 1. Filter significant DEGs (already annotated)
coad_sig <- coad_degs_annot[abs(coad_degs_annot$logFC) > 2 & coad_degs_annot$PValue < 0.01, ]
coad_sig <- coad_sig[!is.na(coad_sig$hgnc_symbol), ]

read_sig <- read_degs_annot[abs(read_degs_annot$logFC) > 2 & read_degs_annot$PValue < 0.01, ]
read_sig <- read_sig[!is.na(read_sig$hgnc_symbol), ]

# 2. Find common genes
common_genes <- intersect(coad_sig$hgnc_symbol, read_sig$hgnc_symbol)
length(common_genes)  # check number of common DEGs

# 3. Subset only common genes and remove duplicates
coad_common <- coad_sig[coad_sig$hgnc_symbol %in% common_genes, ]
coad_common <- coad_common[!duplicated(coad_common$hgnc_symbol), c("hgnc_symbol","logFC","PValue","FDR")]

read_common <- read_sig[read_sig$hgnc_symbol %in% common_genes, ]
read_common <- read_common[!duplicated(read_common$hgnc_symbol), c("hgnc_symbol","logFC","PValue","FDR")]

# 4. Merge COAD and READ info
common_DEGs <- merge(
  coad_common,
  read_common,
  by = "hgnc_symbol",
  suffixes = c("_COAD","_READ")
)
colnames(common_DEGs)[1] <- "GeneSymbol"

# 5. Export merged common DEG table
write.csv(common_DEGs, "TCGA_COAD_READ_common_DEGs44444.csv", row.names = FALSE)

# 6. Quick check
head(common_DEGs)
