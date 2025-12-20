# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
# Differential expression analysis with limma
rm(list=ls())
install.packages("umap")
library(GEOquery)
library(limma)
library(umap)


urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE70073", "file=GSE70073_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
View(tbl)
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID




















*******************************************************
# ONE-CELL SOLUTION for GSE137327 (RNA-seq)

library(GEOquery)
library(R.utils)

# Download supplementary files
getGEOSuppFiles("GSE137327")

# Unzip all files
files <- list.files("GSE137327", full.names = TRUE)
sapply(files, function(f) if (grepl("\\.gz$", f)) gunzip(f, overwrite = TRUE))

# Read expression matrix (adjust filename if needed)
expr_file <- list.files("GSE137327", pattern = "count|FPKM|TPM", full.names = TRUE)[1]
expr <- read.table(expr_file, header = TRUE, row.names = 1, sep = "\t")

# Check expression data
dim(expr)
head(expr)
******************************************


# Load series and platform data from GEO
gset <- getGEO("GSE137327", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL23227", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Make proper column names to match top table
fvarLabels(gset) <- make.names(fvarLabels(gset))

# -------------------------------
# Group membership for all samples
# 9 tumor (T1–T9) and 9 control (C1–C9)
# -------------------------------
# "1" = Colon cancer, "0" = Control
gsms <- "101010101010101010"   # 18 samples: alternating T and C
sml <- strsplit(gsms, split = "")[[1]]

# Log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}
****************
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
LogC <- !any(is.na(qx)) && ((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0))
if (LogC) { ex[ex <= 0 | is.na(ex)] <- NaN; exprs(gset) <- log2(ex) }


# Assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Cancer", "Control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

# Remove missing values
gset <- gset[complete.cases(exprs(gset)), ]

# Fit linear model
fit <- lmFit(gset, design)

# Set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)

# Compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = 250)
colnames(tT)

# View and export results
View(tT)
tT1 <- subset(tT, select = c("ID","adj.P.Val","P.Value","t","B","logFC",
                             "GenBank.Accession","Platform_SPOTID",
                             "Gene.symbol","Gene.title"))
View(tT1)
write.table(tT1, file = stdout(), row.names = FALSE, sep = "\t")

################################################################
# Visualization and QC
tT2 <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
hist(tT2$adj.P.Val, col = "red", border = "white",
     xlab = "P-adjusted value", ylab = "Number of genes",
     main = "P-adj value distribution")

dT <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05, lfc = 0)
vennDiagram(dT, circle.col = palette())

t.good <- which(!is.na(fit2$F))
qqt(fit2$t[t.good], fit2$df.total[t.good], main = "Moderated t statistic")

ct <- 1
volcanoplot(fit2, coef = ct, main = colnames(fit2)[ct], pch = 20,
            highlight = length(which(dT[,ct] != 0)), names = rep('+', nrow(fit2)))

plotMD(fit2, column = ct, status = dT[,ct], legend = FALSE, pch = 20, cex = 1)
abline(h = 0)

################################################################
# General expression data analysis
ex <- exprs(gset)

ord <- order(gs)
palette(c("#1B9E77", "#7570B3"))
par(mar = c(7,4,2,1))
title <- paste("GSE137327", "/", annotation(gset), sep = "")
boxplot(ex[,ord], boxwex = 0.6, notch = TRUE, main = title,
        outline = FALSE, las = 2, col = gs[ord])
legend("topleft", groups, fill = palette(), bty = "n")

par(mar = c(4,4,2,1))
title <- paste("GSE137327", "/", annotation(gset), " value distribution", sep = "")
plotDensities(ex, group = gs, main = title, legend = "topright")

# UMAP plot
ex <- na.omit(ex)
ex <- ex[!duplicated(ex), ]
ump <- umap(t(ex), n_neighbors = 12, random_state = 123)
par(mar = c(3,3,2,6), xpd = TRUE)
plot(ump$layout, main = "UMAP plot, nbrs=12", xlab = "", ylab = "",
     col = gs, pch = 20, cex = 1.5)
legend("topright", inset = c(-0.15,0), legend = levels(gs), pch = 20,
       col = 1:nlevels(gs), title = "Group", pt.cex = 1.5)

install.packages("ggrepel")
library(ggrepel)
library(ggplot2)

df <- data.frame(UMAP1 = ump$layout[,1],
                 UMAP2 = ump$layout[,2],
                 Label = rownames(ump$layout))

ggplot(df, aes(x = UMAP1, y = UMAP2)) +
  geom_point() +
  geom_text_repel(aes(label = Label), size = 3) +
  theme_minimal()

plotSA(fit2, main = "Mean variance trend, GSE137327")