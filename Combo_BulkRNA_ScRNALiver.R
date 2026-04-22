####Load Libraries###
library(Matrix)
library(Seurat)
library(SeuratObject)
library(BiocManager)
library(SingleCellExperiment)
library(scDblFinder)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(fields)
library(presto)
library(SingleR)
library(celldex)
library(dplyr)
###Import files###
###Single Cell RNAseq normal data###
ScRNAseqMatrix_obj <- ReadMtx(mtx = "GSE189175_counts_matrix.mtx.gz", 
                          features ="GSE189175_counts_features.tsv.gz",
                          cells = "GSE189175_counts_barcodes.tsv.gz")
View(ScRNAseqMatrix_obj)
###Create Seurat Object###
Seurat_mtx <- CreateSeuratObject(counts = ScRNAseqMatrix_obj,
                                 project = "HCC_scRNA",
                                 min.cells = 3,
                                 min.features = 200)
###QC and filtration###
View(Seurat_mtx@meta.data)
Seurat_mtx$sample <- rownames(Seurat_mtx@meta.data)
####Split Sample Column###

Seurat_mtx@meta.data <- separate(Seurat_mtx@meta.data, col = 'sample',
                                 into = c('Patient','Barcode'), 
                                 sep = '_')
###Split Barcode Column###
Seurat_mtx@meta.data <- separate(Seurat_mtx@meta.data, col = 'Barcode',
                                 into = c('Tissue','Barcode'), 
                                 sep = ';')

View(Seurat_mtx@meta.data)
####Calculate mitochondrial, ribosomal, and hemoglobin percentages###
Seurat_mtx[["percent.mt"]] <- PercentageFeatureSet(Seurat_mtx, 
                                                   pattern = "^MT-")
View(Seurat_mtx@meta.data)####Are low due to downloaded counts being already processed
x11()
VlnPlot(Seurat_mtx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)
FeatureScatter(Seurat_mtx, feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
Seurat_mtx[["percent.ribo"]] <- PercentageFeatureSet(
  Seurat_mtx,
  pattern = "^RPL|^RPS"
)
View(Seurat_mtx@meta.data)
VlnPlot(Seurat_mtx, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"),
        ncol = 3)

Seurat_mtx[["percent.hb"]] <- PercentageFeatureSet(
  Seurat_mtx,
  pattern = "^HB[AB]"
)
View(Seurat_mtx@meta.data)
VlnPlot(Seurat_mtx, features = c("nFeature_RNA", "nCount_RNA", "percent.hb"),
        ncol = 3)
head(rownames(Seurat_mtx), 20)
View(Seurat_mtx)
#filtering ###
Seurat_mtx_filtered <- subset(Seurat_mtx, subset = nCount_RNA >200 
                              & nFeature_RNA<7000)
dim(Seurat_mtx)
dim(Seurat_mtx_filtered) #how many cells left

###Normalize Data and find features after filtering###
Seurat_mtx_filtered <- NormalizeData(object = Seurat_mtx_filtered)
Seurat_mtx_filtered <- FindVariableFeatures(object = Seurat_mtx_filtered)
top10 <- head(VariableFeatures(Seurat_mtx_filtered), 10)
top10
x11()
VarFeatPlot <- VariableFeaturePlot(Seurat_mtx_filtered)
LabelPoints(plot = VarFeatPlot, points = top10, repel = T)


Seurat_mtx_filtered<- ScaleData(object = Seurat_mtx_filtered)
Seurat_mtx_filtered <- RunPCA(object = Seurat_mtx_filtered, 
                              features = 
                                VariableFeatures(obj=Seurat_mtx_filtered))

print(Seurat_mtx_filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(Seurat_mtx_filtered, dims = 5, cells = 500, balanced = TRUE)

ElbowPlot(Seurat_mtx_filtered) ###Elbow plot made to determine appropriate dimensions###
Seurat_mtx_filtered <- FindNeighbors(object = Seurat_mtx_filtered, dims = 1:15)
Seurat_mtx_filtered <- FindClusters(object = Seurat_mtx_filtered,
                                    resolution = c(0.1,0.3,0.5,0.7,1))
View(Seurat_mtx_filtered@meta.data)
Seurat_mtx_filtered <- RunUMAP(object = Seurat_mtx_filtered, dims = 1:15)
DimPlot(Seurat_mtx_filtered, group.by = "RNA_snn_res.0.1", label = T, 
        reduction = "umap")
DimPlot(Seurat_mtx_filtered, group.by = "RNA_snn_res.0.3", label = T)
DimPlot(Seurat_mtx_filtered, group.by = "RNA_snn_res.0.5", label = T)
DimPlot(Seurat_mtx_filtered, group.by = "RNA_snn_res.0.7", label = T)
DimPlot(Seurat_mtx_filtered, group.by = "RNA_snn_res.1", label = T)

head(Seurat_mtx_filtered@meta.data)

#run scDblFinder from Bioconductor since DoubletFinder Incompatible with Seurat_v5

sce <- as.SingleCellExperiment(Seurat_mtx_filtered)
sce$sampleID <- paste(sce$Patient, sce$Tissue, sep = "_") 
table(sce$sampleID)
sce <- scDblFinder(sce, samples = sce$orig.ident)
table(sce$scDblFinder.class)
###Remove Doublets###
sce <- scDblFinder(sce, samples = sce$sampleID)
table(sce$scDblFinder.class)
sce <- sce[, sce$scDblFinder.class == "singlet"]
table(sce$scDblFinder.class)
Seurat_mtx_filtered_nodb <- as.Seurat(sce)
ncol(Seurat_mtx_filtered_nodb)
names(Seurat_mtx_filtered_nodb@reductions)
# Clear stale reductions
Seurat_mtx_filtered_nodb[["PCA"]] <- NULL
Seurat_mtx_filtered_nodb[["UMAP"]] <- NULL
names(Seurat_mtx_filtered_nodb@reductions)
# Confirm they're gone
names(Seurat_mtx_filtered_nodb@reductions)  # should return character(0)
###Rerun Evalutation on dataset with no doublets###
Seurat_mtx_filtered_nodb <- NormalizeData(Seurat_mtx_filtered_nodb)
Seurat_mtx_filtered_nodb <- FindVariableFeatures(Seurat_mtx_filtered_nodb)
Seurat_mtx_filtered_nodb <- ScaleData(Seurat_mtx_filtered_nodb)
Seurat_mtx_filtered_nodb <- RunPCA(Seurat_mtx_filtered_nodb,
                                   features = VariableFeatures(Seurat_mtx_filtered_nodb))
ElbowPlot(Seurat_mtx_filtered_nodb)
Seurat_mtx_filtered_nodb <- FindNeighbors(Seurat_mtx_filtered_nodb, dims = 1:15)
Seurat_mtx_filtered_nodb <- FindClusters(Seurat_mtx_filtered_nodb,
                                         resolution = c(0.1, 0.3, 0.5, 0.7, 1))

Seurat_mtx_filtered_nodb <- RunUMAP(Seurat_mtx_filtered_nodb, dims = 1:15)
Idents(Seurat_mtx_filtered_nodb) <- "RNA_snn_res.0.5"
###Dimension plots and Heatmaps of newly Evaluated data###
DimPlot(Seurat_mtx_filtered_nodb, label = TRUE, reduction = "umap") +
  xlab("UMAP_1") + ylab("UMAP_2")


DimPlot(Seurat_mtx_filtered_nodb, group.by = "RNA_snn_res.0.1", label = T, 
        reduction = "umap") + xlab("Umap_1") + ylab("Umap_2")

DimPlot(Seurat_mtx_filtered_nodb, group.by = "RNA_snn_res.0.3", label = T, 
        reduction = "umap") + xlab("Umap_1") + ylab("Umap_2")

DimPlot(Seurat_mtx_filtered_nodb, group.by = "RNA_snn_res.0.5", label = T, 
        reduction = "umap") + xlab("Umap_1") + ylab("Umap_2")


DimPlot(Seurat_mtx_filtered_nodb, group.by = "RNA_snn_res.0.7", label = T, 
        reduction = "umap") + xlab("Umap_1") + ylab("Umap_2")

DimPlot(Seurat_mtx_filtered_nodb, group.by = "RNA_snn_res.1", label = T, 
        reduction = "umap") + xlab("Umap_1") + ylab("Umap_2")

DimPlot(Seurat_mtx_filtered_nodb, group.by = "Tissue",reduction = "umap", 
        label = T) + xlab("Umap_1") + ylab("Umap_2")
DimPlot(Seurat_mtx_filtered_nodb, group.by = "Patient",reduction = "umap", 
        label = T) + xlab("Cluster_1") + ylab("Cluster_2")
DoHeatmap(Seurat_mtx_filtered_nodb, features = top10)

class(top10)
str(top10)
ncol(Seurat_mtx_filtered)
ncol(Seurat_mtx_filtered_nodb)
###Find Cell Markers###
markers_pos<- FindAllMarkers(
  Seurat_mtx_filtered_nodb,
  only.pos = T,
  min.pct = 0.1,
  logfc.threshold = 0.25)
head(markers_pos)


levels(Idents(Seurat_mtx_filtered_nodb))
head(Idents(Seurat_mtx_filtered_nodb))
table(sce$scDblFinder.class)  # if sce still exists
names(Seurat_mtx_filtered_nodb@reductions)



top10 <- head(VariableFeatures(Seurat_mtx_filtered_nodb), 10)

# View top markers per cluster
markers_pos %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  View()

####Top Marker Per Cluster###
top10_markers <- markers_pos %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
head(top10_markers)
table(top10_markers$cluster)
x11()
DoHeatmap(Seurat_mtx_filtered_nodb, features = top10_markers$gene)

# Get reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Run annotation
sce_for_singler <- as.SingleCellExperiment(Seurat_mtx_filtered_nodb)
singler_results <- SingleR(
  test = sce_for_singler,
  ref = ref,
  labels = ref$label.main
)

Seurat_mtx_filtered_nodb$singler_labels <- singler_results$labels
x11()
DimPlot(Seurat_mtx_filtered_nodb, 
        group.by = "singler_labels", 
        reduction = "umap",
        label = TRUE,
        repel = TRUE) +
  xlab("UMAP_1") + ylab("UMAP_2")

singler_results$scores
table(singler_results$labels)
x11()
plotScoreHeatmap(singler_results)
plotDeltaDistribution(singler_results)
tab <- table(Assigned=singler_results$labels, Clusters=Seurat_mtx_filtered_nodb$seurat_clusters)
tab
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))
######Bulk RNA Seq Analysis#####
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(remotes)
library(fields)
library(limma)
library(pheatmap)
####Import Bulk RNA Seq Data####
TCGARNAseq <- read.table("TCGA.LIHC.sampleMap-HiSeqV2", 
                         header = T, sep = "\t", row.names = 1)
head(TCGARNAseq)
colnames(TCGARNAseq)
TCGARNAseqmeta <- read.table("TCGA.LIHC.sampleMap-LIHC_clinicalMatrix", 
                             sep = "\t", header = T, fill = T, row.names = 1)
dim(TCGARNAseq)
dim(TCGARNAseqmeta)
###Fixing format to match columns and rows###
colnames(TCGARNAseq) <- gsub("\\.", "-", colnames(TCGARNAseq))
all(colnames(TCGARNAseq) %in% rownames(TCGARNAseqmeta))
TCGARNAseqmeta <- TCGARNAseqmeta[colnames(TCGARNAseq), ]
all(rownames(TCGARNAseqmeta) == colnames(TCGARNAseq))
##Same order of columns and rows in both meta data and sequencing data check##
all(rownames(TCGARNAseqmeta) == colnames(TCGARNAseq))
TCGARNAseqmeta<- TCGARNAseqmeta[colnames(TCGARNAseq), ]
head(TCGARNAseqmeta)
sum(is.na(rownames(TCGARNAseqmeta)))
any(is.na(TCGARNAseqmeta[,1]))
####Remove Recurring Tumor Sample Since there is only 2###
TCGARNAseqmeta_filtered <- TCGARNAseqmeta[
  TCGARNAseqmeta$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"),
]
table(TCGARNAseqmeta_filtered$sample_type)
TCGARNAseq_filtered <- TCGARNAseq[, rownames(TCGARNAseqmeta_filtered)]
dim(TCGARNAseq_filtered)
###Realign columns and rows after removal of Recurring Tumor###
TCGARNAseqmeta_filtered <- TCGARNAseqmeta[
  TCGARNAseqmeta$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"),
]

TCGARNAseq_filtered <- TCGARNAseq[, rownames(TCGARNAseqmeta_filtered)]

TCGARNAseqmeta_filtered <- TCGARNAseqmeta_filtered[colnames(TCGARNAseq_filtered), ]
###Perform Differential Expression###
group <- factor(TCGARNAseqmeta_filtered$sample_type)
group <- relevel(group, ref = "Solid Tissue Normal")
design <- model.matrix(~ group)


fit <- lmFit(TCGARNAseq_filtered, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = 2, number = Inf)
####Verifying Data Transformation####
summary(as.numeric(as.matrix(TCGARNAseq_filtered)))

#####Volcano Plot###
x11()
plot(results$logFC,
     -log10(results$adj.P.Val),
     pch = 20,
     xlab = "Log2 Fold Change",
     ylab = "-Log10 Adjusted P-value",
     main = "Volcano Plot: Tumor vs Normal")
abline(h = -log10(0.05), col = "red")
abline(v = c(-1, 1), col = "blue")
sig_genes <- results[results$adj.P.Val < 0.05, ]
####PCA plot###
pca <- prcomp(t(TCGARNAseq_filtered))

plot(pca$x[,1], pca$x[,2],
     col = as.factor(TCGARNAseqmeta_filtered$sample_type),
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of TCGA LIHC Samples")
top_genes <- rownames(results)[1:50]
x11()
pheatmap(TCGARNAseq_filtered[top_genes, ])
#####Integration####
#####################
#####Load Libraries####
library(dplyr)
library(scrapper)
library(celldex)
library(SingleR)
library(VennDiagram)
### Identify bulk DEGs
deg_bulk <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]

bulk_sig  <- rownames(deg_bulk)
bulk_up   <- rownames(results[results$logFC > 1 & results$adj.P.Val < 0.05, ])
bulk_down <- rownames(results[results$logFC < -1 & results$adj.P.Val < 0.05, ])

head(rownames(deg_bulk))
head(markers_pos$gene)
length(markers_pos$gene)
length(unique(markers_pos$gene))
length(rownames(deg_bulk))
### Find overlap between bulk DEGs and scRNA markers
overlap_genes <- intersect(rownames(deg_bulk), markers_pos$gene)
length(overlap_genes)
overlap_markers <- markers_pos[
  markers_pos$gene %in% overlap_genes,
]
table(overlap_markers$cluster)

length(overlap_genes)
table(overlap_markers$cluster)
### Cluster enrichment###
sort(table(overlap_markers$cluster), decreasing = TRUE)
top_cluster <- names(sort(table(overlap_markers$cluster), decreasing = TRUE))[1]

top_cluster_genes <- overlap_markers[
  overlap_markers$cluster == top_cluster,
]
table(overlap_markers$cluster)
colnames(markers_pos)
head(top_cluster_genes[order(-top_cluster_genes$avg_log2FC), ])
colnames(markers_pos)

Seurat_mtx_filtered_nodb$singler_labels
Idents(Seurat_mtx_filtered_nodb)
as.data.frame(table(
  Cluster = Idents(Seurat_mtx_filtered_nodb),
  CellType = Seurat_mtx_filtered_nodb$singler_labels
))
subset(
  as.data.frame(
    table(
      Cluster = Idents(Seurat_mtx_filtered_nodb),
      CellType = Seurat_mtx_filtered_nodb$singler_labels
    )
  ),
  Cluster == 6
)
head(top_cluster_genes[order(-top_cluster_genes$avg_log2FC), ])

up_overlap <- intersect(bulk_up, markers_pos$gene)
down_overlap <- intersect(bulk_down, markers_pos$gene)

###Plot intersections###
length(rownames(deg_bulk))              # bulk DEGs
length(unique(markers_pos$gene))        # scRNA markers (unique!)
length(intersect(rownames(deg_bulk), markers_pos$gene))  # overlap
bulk_only <- setdiff(rownames(deg_bulk), markers_pos$gene)
scrna_only <- setdiff(unique(markers_pos$gene), rownames(deg_bulk))
overlap <- intersect(rownames(deg_bulk), markers_pos$gene)
length(bulk_only) + length(overlap) == length(rownames(deg_bulk))
length(scrna_only) + length(overlap) == length(unique(markers_pos$gene))
length(bulk_only)
length(scrna_only)
length(overlap)
x11()
venn.plot <- venn.diagram(
  x = list(
    Bulk_DEGs = rownames(deg_bulk),
    scRNA_markers = unique(markers_pos$gene)
  ),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)


grid::grid.draw(venn.plot)

cluster_counts <- table(overlap_markers$cluster)

barplot(
  sort(cluster_counts, decreasing = TRUE),
  las = 2,
  col = "steelblue",
  main = "Overlap Genes per Cluster",
  ylab = "Number of overlapping genes"
)


df <- as.data.frame(table(overlap_markers$cluster))

ggplot(df, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster", y = "Overlap Gene Count",
       title = "Bulk–scRNA Overlap by Cluster") +
  theme_minimal()

FeaturePlot(
  Seurat_mtx_filtered_nodb,
  features = c("HAMP", "FAM99A", "TMPRSS15"),
  reduction = "umap"
)

Seurat_mtx_filtered_nodb <- AddModuleScore(
  Seurat_mtx_filtered_nodb,
  features = list(overlap_genes),
  name = "OverlapScore"
)

FeaturePlot(
  Seurat_mtx_filtered_nodb,
  features = "OverlapScore1"
)

barplot(
  sort(table(overlap_markers$cluster), decreasing = TRUE),
  las = 2,
  main = "Overlap Genes per Cluster"
)





