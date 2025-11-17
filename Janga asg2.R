# Assignment 2 – SARS-CoV-2 miRNA Differential Expression

# Load required packages
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)

setwd("C:/Users/charu/OneDrive/Desktop/outputs")
# 1. Load Counts Matrix

counts <- read.csv(
  "C:/Users/charu/OneDrive/Desktop/counts final csv.csv",
  row.names = 1,
  check.names = FALSE
)

dim(counts)       

# 2. Metadata 

meta <- data.frame(
  condition = factor(c(
    "SARS","SARS","SARS",
    "Mock","Mock","Mock",
    "SARS","SARS","SARS",
    "Mock","Mock","Mock"
  )),
  
  time = factor(c(
    "72h","72h","72h",
    "72h","72h","72h",
    "24h","24h","24h",
    "24h","24h","24h"
  )),
  
  row.names = c(
    "SRR22269872",
    "SRR22269873",
    "SRR22269874",
    "SRR22269875",
    "SRR22269876",
    "SRR22269877",
    "SRR22269878",
    "SRR22269879",
    "SRR22269880",
    "SRR22269881",
    "SRR22269882",
    "SRR22269883"
  )
)

# Make sure metadata matches counts
all(rownames(meta) == colnames(counts))   

# 3. Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition + time
)


# 4. Run DESeq2
dds <- DESeq(dds)

# 5. PCA 
vsd <- vst(dds)
plotPCA(vsd, intgroup = c("condition","time"))


# 6. Differential Expression Analyses

### 6.1 SARS vs Mock (all timepoints)
res_sars_vs_mock <- results(dds, contrast = c("condition","SARS","Mock"))

# Shrink LFC using normal shrinkage
res_sars_vs_mock <- lfcShrink(
  dds,
  contrast = c("condition","SARS","Mock"),
  res = res_sars_vs_mock,
  type = "normal"
)

DEGs_sars_vs_mock <- res_sars_vs_mock[
  which(res_sars_vs_mock$padj < 0.05 &
          abs(res_sars_vs_mock$log2FoldChange) > 1),
]
# Save results
write.csv(as.data.frame(res_sars_vs_mock),
          "DEG_SARS_vs_Mock_full.csv")
write.csv(as.data.frame(DEGs_sars_vs_mock),
          "DEG_SARS_vs_Mock_significant.csv")

### 6.2 SARS 24h vs SARS 72h

dds_sars <- dds[, dds$condition == "SARS"]
dds_sars$time <- factor(dds_sars$time, levels = c("24h","72h"))
design(dds_sars) <- ~ time
dds_sars <- DESeq(dds_sars)

res_sars_24_vs_72 <- results(dds_sars, contrast = c("time","24h","72h"))

DEGs_sars_24_vs_72 <- res_sars_24_vs_72[
  which(res_sars_24_vs_72$padj < 0.05 &
          abs(res_sars_24_vs_72$log2FoldChange) > 1),
]

write.csv(as.data.frame(res_sars_24_vs_72),
          "DEG_SARS_24h_vs_72h_full.csv")
write.csv(as.data.frame(DEGs_sars_24_vs_72),
          "DEG_SARS_24h_vs_72h_significant.csv")


# 7. GO Enrichment Analysis

### Function to convert Ensembl to Entrez safely

convert_to_entrez <- function(gene_list){
  gene_list <- sub("\\..*", "", gene_list)   
  mapIds(org.Hs.eg.db,
         keys = gene_list,
         column = "ENTREZID",
         keytype = "ENSEMBL",
         multiVals = "first") %>% na.omit()
}

### 7.1 GO Analysis - SARS vs Mock
entrez1 <- convert_to_entrez(rownames(DEGs_sars_vs_mock))

entrez1_valid <- bitr(entrez1,
                      fromType="ENTREZID",
                      toType="SYMBOL",
                      OrgDb=org.Hs.eg.db)$ENTREZID

go_sars_vs_mock <- enrichGO(
  gene         = entrez1_valid,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

write.csv(as.data.frame(go_sars_vs_mock),
          "GO_SARS_vs_Mock.csv", row.names = FALSE)

### 7.2 GO Analysis - SARS 24h vs 72h
entrez2 <- convert_to_entrez(rownames(DEGs_sars_24_vs_72))

entrez2_valid <- bitr(entrez2,
                      fromType="ENTREZID",
                      toType="SYMBOL",
                      OrgDb=org.Hs.eg.db)$ENTREZID

go_sars_24_72 <- enrichGO(
  gene         = entrez2_valid,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

write.csv(as.data.frame(go_sars_24_72),
          "GO_SARS_24h_vs_72h.csv", row.names = FALSE)



# 8. Save Normalized Expression Matrix

rld <- rlog(dds)
norm_mat <- assay(rld)
write.csv(norm_mat, "Normalized_miRNA_expression_matrix.csv")


# DEG Summary

library(dplyr)

setwd("C:/Users/charu/OneDrive/Desktop/outputs")


# 1. Load your DEG result files

deg_mock_sars <- read.csv("DEG_SARS_vs_Mock_significant.csv", row.names = 1)
deg_24_72     <- read.csv("DEG_SARS_24h_vs_72h_significant.csv", row.names = 1)

# Convert to data frame (already is, but ensures consistency)
deg_mock_sars <- as.data.frame(deg_mock_sars)
deg_24_72     <- as.data.frame(deg_24_72)

# 2. Summarizing DEGs

summarize_degs <- function(df){
  data.frame(
    Total_DEGs    = nrow(df),
    Upregulated   = sum(df$log2FoldChange > 0, na.rm = TRUE),
    Downregulated = sum(df$log2FoldChange < 0, na.rm = TRUE),
    Mean_log2FC   = round(mean(df$log2FoldChange, na.rm = TRUE), 3),
    Median_log2FC = round(median(df$log2FoldChange, na.rm = TRUE), 3),
    Min_padj      = signif(min(df$padj, na.rm = TRUE), 3),
    Max_padj      = signif(max(df$padj, na.rm = TRUE), 3)
  )
}

summary_mock_sars <- summarize_degs(deg_mock_sars)
summary_sars_24_72 <- summarize_degs(deg_24_72)

print("Mock vs SARS-CoV-2 Summary:")
print(summary_mock_sars)

print("SARS-CoV-2 24h vs 72h Summary:")
print(summary_sars_24_72)

# Save summary tables
write.csv(summary_mock_sars, "Summary_DEGs_Mock_vs_SARS.csv")
write.csv(summary_sars_24_72, "Summary_DEGs_SARS_24h_vs_72h.csv")

# 3. Extract Up/Down Regulated Genes

# Mock vs SARS
up_mock_sars <- rownames(deg_mock_sars[deg_mock_sars$log2FoldChange > 0, ])
down_mock_sars <- rownames(deg_mock_sars[deg_mock_sars$log2FoldChange < 0, ])

# 24h vs 72h
up_24_72 <- rownames(deg_24_72[deg_24_72$log2FoldChange > 0, ])
down_24_72 <- rownames(deg_24_72[deg_24_72$log2FoldChange < 0, ])

# Print them
up_mock_sars
down_mock_sars
up_24_72
down_24_72

# Save gene lists
write.csv(up_mock_sars, "Upregulated_Mock_vs_SARS.csv", row.names = FALSE)
write.csv(down_mock_sars, "Downregulated_Mock_vs_SARS.csv", row.names = FALSE)

write.csv(up_24_72, "Upregulated_SARS_24h_vs_72h.csv", row.names = FALSE)
write.csv(down_24_72, "Downregulated_SARS_24h_vs_72h.csv", row.names = FALSE)
