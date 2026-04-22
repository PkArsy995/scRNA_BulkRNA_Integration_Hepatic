# Integration of Single Cell RNA Sequencing and Bulk RNA Sequencing Data in the Context of Liver Cancer and Normal Liver Cells
## Identifies Common Genes Between Normal functioning Liver Cells and Cancerous Liver Cells

Using public datasets from 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189175 and https://xenabrowser.net/datapages/?cohort=TCGA%20Liver%20Cancer%20(LIHC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 

single-cell RNA-sequencing data were integrated with bulk RNA-seq from a Liver Cells in a normal vs. disease context. The aim is to investigate how bulk transcriptional signals relate to underlying cellular composition or cell-type–specific transcriptional programs in terms of the specific genes that overlap between the cancerous liver cells and normal functioning liver cells.

###Visual Samples of Steps included in Script###
###QC filtration###
<img width="189" height="183" alt="image" src="https://github.com/user-attachments/assets/4db2257c-f77a-4589-966d-7e41e767ea8c" />
Volcano plot of Percentage of mitochondria in cell

<img width="217" height="158" alt="image" src="https://github.com/user-attachments/assets/47cab734-f449-447f-b050-b25748209ba4" />
Feature scatter plot to determine the quality of cells

<img width="261" height="190" alt="image" src="https://github.com/user-attachments/assets/6c1cf570-4ad5-4019-b3f8-73c2ed5bf476" />
Outlier genes expressed in a Variable Feature plot

<img width="214" height="147" alt="image" src="https://github.com/user-attachments/assets/7cf54353-9b11-4daa-9ce2-5dfee5c95640" />
Elbow plot to determine principle components that should be kept for downstream analysis before doublet removal

<img width="314" height="292" alt="image" src="https://github.com/user-attachments/assets/2b35fac4-bcce-4ff0-bacd-4061e0b5e29d" />
Dimension plot before doublet removal

<img width="233" height="232" alt="image" src="https://github.com/user-attachments/assets/2c043c9c-942b-49d0-b4fa-9d35a42a03bd" />
Dimension plot after doublet removal

<img width="261" height="242" alt="image" src="https://github.com/user-attachments/assets/45c1a2e8-e42b-4d2f-8e1b-42b4493750b4" />
Dimension plot after doublet removal labelled with cell types

<img width="255" height="205" alt="image" src="https://github.com/user-attachments/assets/3562288c-5067-46dc-8d1d-11c0b3f157eb" />
Volcano plot of Tumor vs. Normal Genes in Bulk RNA Sequencing

<img width="321" height="184" alt="image" src="https://github.com/user-attachments/assets/0b3cc1e5-4729-49d6-b18c-82ba92eac4ab" />
PCA plot of Tumor cells (black) and normal cells (Red)

<img width="468" height="298" alt="image" src="https://github.com/user-attachments/assets/6a11060e-565f-4124-9ef1-23e0f9768a4a" />
Overlap of ScRNA and Bulk RNAseq Data Bar graph

<img width="468" height="298" alt="image" src="https://github.com/user-attachments/assets/a48bf8c9-5c58-4d43-b346-89c9f2501ecc" />

Examples of Genes overlap expressed in a dimension plot





###Instructions###
1) Install latest versions of R and RStudio depending on operating system used (Mac OS or Windows) via: https://cran.r-project.org
2) Open RStudio and install or initiate the following packages from Packages tab at bottom right hand corner:
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
   library(remotes)
   library(limma)
   library(pheatmap)
   library(scrapper)
   library(VennDiagram)
3) Set working directory and run above RScript




