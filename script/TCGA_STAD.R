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
stad <- get_tcga("TCGA-STAD")

# Expression matrices
stad_counts <- assay(stad)

# Summary table
# Use `sample_type` which is more consistent across TCGA projects
stad_type <- colData(stad)$shortLetterCode
if(is.null(stad_type)) stad_type <- colData(stad)$sample_type

data.frame(
  Dataset = "TCGA-STAD",
  No_of_Normal = sum(stad_type %in% c("NT", "Solid Tissue Normal")),
  No_of_Tumor  = sum(stad_type %in% c("TP", "Primary Tumor")),
  PlatformID  = if("platform" %in% colnames(colData(stad))) unique(colData(stad)$platform)[1] else NA,
  No_of_Row   = nrow(stad_counts)
)

# Keep genes with counts per million > 10 in at least, say, 10% of samples
keep_stad <- rowSums(cpm(stad_counts) > 10) >= 0.1*ncol(stad_counts)
stad_counts <- stad_counts[keep_stad, ]
View(stad_counts)

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
  degs$Significant <- abs(degs$logFC) > 1 & degs$PValue < 0.01
  return(degs)
}

stad_degs <- dge_analysis(stad_counts, stad_type)

sum(stad_degs$Significant)

# Separate Up- and Down-regulated genes
#STAD
stad_up   <- stad_degs[stad_degs$logFC >  1 & stad_degs$PValue < 0.01, ]
stad_down <- stad_degs[stad_degs$logFC < -1 & stad_degs$PValue < 0.01, ]
nrow(stad_down)
nrow(stad_up)
nrow(stad_degs)

### Enhanced Volcano
EnhancedVolcano(stad_degs,
                lab = rownames(stad_degs),
                x = "logFC",
                y = "PValue",
                FCcutoff = 1,
                pCutoff = 0.01,
                title = "TCGA-STAD Volcano Plot")


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
stad_degs_annot <- annotate_deg(stad_degs)
stad_up_annot   <- annotate_deg(stad_up)
stad_down_annot <- annotate_deg(stad_down)

sum(is.na(stad_degs_annot$GeneSymbol))


write.csv(stad_degs, "F:\\DATA_RNASeq\\Cancer-Biology\\data\\TCGA-STAD\\TCGA_STAD_all_DEGs.csv")
write.csv(stad_up,   "F:\\DATA_RNASeq\\Cancer-Biology\\data\\TCGA-STAD\\TCGA_STAD_upregulated.csv")
write.csv(stad_down, "F:\\DATA_RNASeq\\Cancer-Biology\\data\\TCGA-STAD\\TCGA_STAD_downregulated.csv")

# 4. Export CSV files
write.csv(stad_degs_annot, "F:\\DATA_RNASeq\\Cancer-Biology\\data\\TCGA-STAD\\TCGA_stad_all_DEGs_gene.csv", row.names = FALSE)
write.csv(stad_up_annot,   "F:\\DATA_RNASeq\\Cancer-Biology\\data\\TCGA-STAD\\TCGA_stad_upregulated_gene.csv", row.names = FALSE)
write.csv(stad_down_annot, "F:\\DATA_RNASeq\\Cancer-Biology\\data\\TCGA-STAD\\TCGA_stad_downregulated_gene.csv", row.names = FALSE)

top100_deg <- stad_degs_annot[order(stad_degs_annot$PValue), ][1:100, ]

top100_up <- stad_degs_annot[
  stad_degs_annot$logFC > 1,
][order(stad_degs_annot$PValue[stad_degs_annot$logFC > 1]), ][1:100, ]

write.csv(top100_deg,  "F:\\DATA_RNASeq\\Cancer-Biology\\data\\TCGA-STAD\\Top100_DEGs.csv")
write.csv(top100_up,   "F:\\DATA_RNASeq\\Cancer-Biology\\data\\TCGA-STAD\\Top100_Upregulated.csv")


