#Libraries and seeds
library(devtools)
library(R.utils)
library(utils)
library(ggplot2)
library(dplyr)
library(patchwork)
library(SeuratObject)
library(Seurat)
library(hdf5r)
library(Signac)
library(BiocManager)
library(GenomeInfoDb)
library(biovizBase)
library(readr)
library(magrittr)
library(Matrix)
library(tidyr)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(limma)
library(lares)
library(vioplot)
library(biomaRt)
library(ape)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(Connectome)
library(gridGraphics)
library(networkD3)
library(proxyC)
library(ShortRead)
library(dendextend)
library(ggthemes)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(CellChat)
library(WGCNA)
library(spacexr)
library(Matrix)
library(data.table)
library(slingshot)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)
library(monocle3)
library(SeuratWrappers)
library(SCORPIUS)
library(goseq)
library(clusterProfiler)
library(enrichplot)
library(RaceID)

library(hdWGCNA)
library(pbacBarcode)
set.seed(123)
print("done")

#########Importing Data##########
#Mixed.csv is the Mixed population and IUE_Data_Rep2
#Replicate1.csv is the Rep1_F population and IUE_Data_Rep1
#Replicate2.csv is the Rep2_M population and IUE_Data_Rep3
IUE_Data<-read.csv("Replicate1.csv", header=TRUE, comment.char = "#", row.names = 1)
IUE_Data<-as.data.frame(t(IUE_Data))
P4_Rep1 <- CreateSeuratObject(counts = IUE_Data, project = "IUE_Data_Rep1", min.cells = 3, min.features = 200)

IUE_Data_2<-read.csv("Mixed.csv", header=TRUE, comment.char = "#", row.names = 1)
IUE_Data_2<-as.data.frame(t(IUE_Data_2))
P4_Orig <- CreateSeuratObject(counts = IUE_Data_2, project = "IUE_Data_P4", min.cells = 2, min.features = 200)

IUE_Data_3<-read.csv("Replicate2.csv", header=TRUE, comment.char = "#", row.names = 1)
IUE_Data_3<-as.data.frame(t(IUE_Data_3))
P4_Rep2 <- CreateSeuratObject(counts = IUE_Data_3, project = "IUE_Data_Rep3", min.cells = 3, min.features = 200)

P4_Obj <- merge(P4_Rep1, y = c(P4_Orig, P4_Rep2), add.cell.ids = c("Rep1_F", "Mixed", "Rep2_M"), project = "Merged_EGFP")

##Start analysis##
P4_Obj <- PercentageFeatureSet(P4_Obj, pattern = "^mt.", col.name = "percent.mt")
VlnPlot(P4_Obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

#visualize features
P4_plot1<-FeatureScatter(P4_Obj, feature1 = "nCount_RNA", feature2 ="percent.mt")
P4_plot2<-FeatureScatter(P4_Obj, feature1 = "nCount_RNA", feature2 ="nFeature_RNA")
P4_plot1 + P4_plot2

#Subset object
P4_Obj <- subset(P4_Obj, subset=nFeature_RNA > 200 & nCount_RNA < 40000 & percent.mt < 25)

P4_Obj <- NormalizeData(P4_Obj)
P4_Obj <- FindVariableFeatures(P4_Obj, selection.method = "vst")
P4_Obj <- ScaleData(P4_Obj, features = rownames(P4_Obj))

#Genes for cell cycle scoring
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019 $g2m.genes

convertHumanGeneList <- function(x){
  require("biomaRt")
  #human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
s.genes <- convertHumanGeneList(s.genes)
g2m.genes <- convertHumanGeneList(g2m.genes)

#Cell cycle scoring
P4_Obj <- CellCycleScoring(P4_Obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(P4_Obj, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)

#SCTransform
P4_Obj <- SCTransform(P4_Obj, method = "glmGamPoi", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE)

P4_top10 <- head(VariableFeatures(P4_Obj), 10)
P4_top_plot1 <- VariableFeaturePlot(P4_Obj)
P4_top_plot2 <- LabelPoints(plot = P4_top_plot1, points = P4_top10, repel=TRUE)
P4_top_plot1 + P4_top_plot2

P4_Obj <- RunPCA(P4_Obj, verbose = FALSE)

ElbowPlot(P4_Obj, ndims=50)

#Generating UMAP
P4_Obj <- RunUMAP(P4_Obj, dims = 1:35, verbose = FALSE)
P4_Obj <- FindNeighbors(P4_Obj, dims = 1:35, verbose = FALSE)
P4_Obj <- FindClusters(P4_Obj, verbose = FALSE)
DimPlot(P4_Obj, label = TRUE) + NoLegend()
DimPlot(P4_Obj, reduction = "umap", label = TRUE, group.by = "orig.ident", pt.size = 0.5) + NoLegend()
DimPlot(P4_Obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 0.5) + NoLegend()

#Identifying top genes of each cluster
cluster0.markers <- FindMarkers((P4_Obj), ident.1=0, min.pct = 0.25)
Cluster0.markers <- row.names(cluster0.markers)

cluster1.markers <- FindMarkers((P4_Obj), ident.1=1, min.pct = 0.25)
Cluster1.markers <- row.names(cluster1.markers)

cluster2.markers <- FindMarkers((P4_Obj), ident.1=2, min.pct = 0.25)
Cluster2.markers <- row.names(cluster2.markers)

cluster3.markers <- FindMarkers((P4_Obj), ident.1=3, min.pct = 0.25)
Cluster3.markers <- row.names(cluster3.markers)

cluster4.markers <- FindMarkers((P4_Obj), ident.1=4, min.pct = 0.25)
Cluster4.markers <- row.names(cluster4.markers)

cluster5.markers <- FindMarkers((P4_Obj), ident.1=5, min.pct = 0.25)
Cluster5.markers <- row.names(cluster5.markers)

cluster6.markers <- FindMarkers((P4_Obj), ident.1=6, min.pct = 0.25)
Cluster6.markers <- row.names(cluster6.markers)

cluster7.markers <- FindMarkers((P4_Obj), ident.1=7, min.pct = 0.25)
Cluster7.markers <- row.names(cluster7.markers)

cluster8.markers <- FindMarkers((P4_Obj), ident.1=8, min.pct = 0.25)
Cluster8.markers <- row.names(cluster8.markers)

cluster9.markers <- FindMarkers((P4_Obj), ident.1=9, min.pct = 0.25)
Cluster9.markers <- row.names(cluster9.markers)

cluster10.markers <- FindMarkers((P4_Obj), ident.1=10, min.pct = 0.25)
Cluster10.markers <- row.names(cluster10.markers)

cluster11.markers <- FindMarkers((P4_Obj), ident.1=11, min.pct = 0.25)
Cluster11.markers <- row.names(cluster11.markers)

cluster12.markers <- FindMarkers((P4_Obj), ident.1=12, min.pct = 0.25)
Cluster12.markers <- row.names(cluster12.markers)

cluster13.markers <- FindMarkers((P4_Obj), ident.1=13, min.pct = 0.25)
Cluster13.markers <- row.names(cluster13.markers)

cluster14.markers <- FindMarkers((P4_Obj), ident.1=14, min.pct = 0.25)
Cluster14.markers <- row.names(cluster14.markers)

cluster15.markers <- FindMarkers((P4_Obj), ident.1=15, min.pct = 0.25)
Cluster15.markers <- row.names(cluster15.markers)

cluster16.markers <- FindMarkers((P4_Obj), ident.1=16, min.pct = 0.25)
Cluster16.markers <- row.names(cluster16.markers)

cluster17.markers <- FindMarkers((P4_Obj), ident.1=17, min.pct = 0.25)
Cluster17.markers <- row.names(cluster17.markers)

cluster18.markers <- FindMarkers((P4_Obj), ident.1=18, min.pct = 0.25)
Cluster18.markers <- row.names(cluster18.markers)

cluster19.markers <- FindMarkers((P4_Obj), ident.1=19, min.pct = 0.25)
Cluster19.markers <- row.names(cluster19.markers)

cluster20.markers <- FindMarkers((P4_Obj), ident.1=20, min.pct = 0.25)
Cluster20.markers <- row.names(cluster20.markers)

#These will print you a jpeg of the top 12 genes for each cluster
jpeg("P4_Obj_Cluster0.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster0.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster1.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster1.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster2.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster2.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster3.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster3.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster4.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster4.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster5.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster5.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster6.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster6.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster7.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster7.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster8.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster8.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster9.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster9.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster10.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster10.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster11.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster11.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster12.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster12.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster13.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster13.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster14.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster14.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster15.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster15.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster16.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster16.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster17.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster17.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster18.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster18.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster19.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster19.markers, 12)))
dev.off()

jpeg("P4_Obj_Cluster20.jpg", width = 1920, height = 1080)
FeaturePlot(P4_Obj, c(head(Cluster20.markers, 12)))
dev.off()

#Highlights cells within each cluster
jpeg("P4_Obj_Cluster0_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 0, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster1_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 1, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster2_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 2, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster3_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 3, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster4_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 4, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster5_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 5, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster6_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 6, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster7_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 7, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster8_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 8, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster9_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 9, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster10_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 10, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster11_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 11, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster12_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 12, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster13_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 13, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster14_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 14, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster15_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 15, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Objj_Cluster16_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 16, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster17_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 17, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster18_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 18, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster19_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 19, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

jpeg("P4_Obj_Cluster20_Highlight.jpg", width = 1920, height = 1080)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5, order = 20, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "red")) + NoLegend()
dev.off()

###   Naming each cluster with their cell-type, these are best guesses based upon markers from earlier plots, 4,6,8 are not necessarily committed but their closest matching "mature" populations are interneurons, oligodendrocytes and astrocytes respectively.
CID0 <- "0 Radial Glia/Adult Neural Stem Cells"
CID1 <- "1 Astrocytes" 
CID2 <- "2 Oligodendrocytes"
CID3 <- "3 Maturing Olfactory Interneurons"
CID4 <- "4 Olfactory Interneuron Committed TIPs"
CID5 <- "5 Mature Olfactory Interneurons"
CID6 <- "6 Oligodendrocyte Committed TIPs"
CID7 <- "7 Excitatory Neurons" 
CID8 <- "8 Tripotent Intermediate Progenitors (TIPs)"
CID9 <- "9 Ependymal-like Radial Glia"
CID10 <- "10 Oligodendrocyte Progenitor Cells (OPCs)"
CID11 <- "11 Radial Glia/Adult Neural Stem Cells"
CID12 <- "12 Maturing Excitatory Neurons"
CID13 <- "13 Astrocyte Subpopulation"
CID14 <- "14 Mature Oligodendrocytes"
CID15 <- "15 Unknown"
CID16 <- "16 Unknown"
CID17 <- "17 Ependymal Cells"
CID18 <- "18 Unknown"
CID19 <- "19 Olfactory Interneuron Subpopulation"
CID20 <- "20 Unknown"

new.cluster.ids <- c(CID0, CID1, CID2, CID3, CID4, CID5, CID6, CID7, CID8, CID9, CID10, CID11, CID12, CID13, CID14, CID15, CID16, CID17, CID18, CID19, CID20)
names(new.cluster.ids) <- levels(P4_Obj)
P4_Obj <- RenameIdents(P4_Obj, new.cluster.ids)
DimPlot(P4_Obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(P4_Obj, reduction = "umap", label = FALSE, pt.size = 0.5, split.by = "orig.ident") + NoLegend()

# Subsetting by useful clusters/removing the unknown clusters
P4_Clean <- subset(P4_Obj, subset = seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 3 | seurat_clusters == 4 | seurat_clusters == 5 | seurat_clusters == 6 | seurat_clusters == 7 | seurat_clusters == 8 | seurat_clusters == 9 | seurat_clusters == 10 | seurat_clusters == 11 | seurat_clusters == 12 | seurat_clusters == 13 | seurat_clusters == 14 | seurat_clusters == 17 | seurat_clusters == 19) 


############################################################
#   #####   #          #      ##      #   ########   ###   #
#  #        #         # #     # #     #   #         #      #
# #         #        #   #    #  #    #   #        #       #
# #         #        #   #    #   #   #   ######    #####  #
# #         #        #   #    #    #  #   #             #  #
#  #        #         # #     #     # #   #            #   #
#   #####   #######    #      #      ##   ########  ###    #
############################################################

#This is the pbacBarcode package used to assign clones. Here we demonstrate its use for calculating the clones in the Mixed population.
#Files:
#pbCellTag1_S1_L001_R1_001.fastq
#pbCellTag_and_Cell_filtered_IUE_Apr8.bam
#These are available on request as they are the filtered outputs of BD's pipeline
# The following is an example of how clones were generated

# To begin testing the package, please install the pbCellTag package using this command in RStudio.  It takes some time to install because of the sample data files.
install.packages("pbacBarcode_1.0.0.tar.gz", repos = NULL, type="source")

# After it has installed, simply load it:
library(pbacBarcode)
# The script below follows the code from the CellTagR vignette (https://github.com/morris-lab/CellTagR)

### Assessment of CellTag Library Complexity via Sequencing ----
# -- 1. Read in the FASTQ sequencing data and extract the CellTags --

# Read in the data file. Note: add your path to the FASTQ file
rhap_path <- "/example_path/pbCellTag1_S1_L001_R1_001.fastq"
rhapsody_obj <- CellTagObject(object_name = "rhapsody_fastq",
                              file_directory = rhap_path,
                              technique = "Rhapsody")

# Extract the CellTag patterns using the CellTagPatternCalling() function. Submit your CellTag version and the scRNAseq technique. You will get a list of the full regex CellTag string, the prefix nucleotides and the suffix nucleotides that are found in the beginning and end of the full CellTag string.
patterns <- CellTagPatternCalling("v3", "rhapsody")
rhapsody_obj <- CTExtraction(pbacBarcode_obj = rhapsody_obj, celltag_version = "v3", patterns)

# -- 2. Count the CellTags and sort based on the occurrence of each CellTag --
rhapsody_obj <- AddCellTagFreqSort(rhapsody_obj)

# Check the stats
head(rhapsody_obj@celltag_freq_stats$v3)

# -- 3. Generation of a whitelist for the CellTag library --
# This code will generate the 99th percentile whitelist used in the paper
rhapsody_obj <- CTWhitelistFiltering(rhapsody_obj,
                                     output_dir = "./rhapsody_whitelist_Mixed_99.csv",
                                     percentile = 0.99,
                                     save_file = T)
# The generated whitelist for each library will be used downstream to filter and clean the single-cell CellTag UMI matrices.

# Now, you will be working with the BAM file and will use that whitelist file for filtering.

### Single-Cell CellTag Extraction and Quantification ----
# In this section, an alternative approach is demonstrated that conducts CellTag extraction, quantification, and generation of UMI count matrices.
# -- 1. Create a CellTag Object --
rm(rhapsody_obj) # this was used to create the whitelist file.
rhap_path <- "/example_path/pbCellTag_and_Cell_filtered_IUE_Apr8.bam"
rhapsody_obj <- CellTagObject("rhapsody_bam", rhap_path, technique = "Rhapsody")

# -- 2. Extract the CellTags from the BAM file --
# Get the CellTag patterns using CellTagPatternCalling(). Extract the dataframe of CellTags, barcodes, and UMIs using the CTExtraction() function. The run-time depends on the size of the BAM file. It needs to read the BAM file (which is the slowest processed). The preceeding steps are quick.
patterns <- CellTagPatternCalling("v3", "rhapsody")
rhapsody_obj <- CTExtraction(pbacBarcode_obj = rhapsody_obj, "v3", celltag_patterns = patterns)

# -- 3. Quantify the CellTag UMI Counts and Generate UMI Count Matrices --
# Generate the sparse count matrix. The generated CellTag UMI count matrices can then be used in for clone identification. In the next section, different data is used to demonstrate how to retrieve clones from CellTag UMI count matrices. You can follow the steps in the next section to retrieve the clones from this dataset.
fpath <- "/example_path/IUE_Cell_Barcodes.tsv"
rhapsody_obj <- CTMatrixCount(rhapsody_obj, fpath)
# The generated CellTag UMI count matrices can then be used in the following steps for clone identification.

### Single-cell CellTag UMI Count Matrix Processing ----
# --- I. Prepare for the data to be collapsed ---
# In this section, starcode will be used to cluster the cells. First, prepare the data to the format that is accepted by starcode. This function accepts two inputs including the CellTag object with raw count matrix generated and a path to where to save the output text file. The output will be a text file with each line containing one sequence to collapse with others. In this function, we concatenate the CellTag with cell barcode and use the combined sequences as input to execute Starcode. The file to be used for Starcode will be stored under the provided directory.
# Generating the collapsing file
rhapsody_obj <- CTForCollapsing(rhapsody_obj, output_file = "/example_path/collapsing_input_Mixed.txt")

# --- II. Run Starcode to cluster CellTags ---
# Following the instruction for Starcode, run the following command to generate the result from starcode. Make sure to change the directories wherever starcode is and where you want to save the resulting file.
# ./starcode -s --print-clusters collapsing_input_Mixed.txt > collapsing_result_Mixed.txt
rhapsody_obj <- CTPostCollapsing(rhapsody_obj, "/example_path/collapsing_result_Rep1_OCT112024.txt")

# -- 2. Binarize the single-cell CellTag UMI count matrix --
# Calling binarization
rhapsody_obj <- MatrixBinarization(rhapsody_obj, 2)
mtx_test <- pbacBarcode::GetMatrixByVersion(rhapsody_obj, "binary_mtx")
celltags_per_cell_whitelisted <- rowSums(mtx_test)
celltags_frequency_whitelisted <- colSums(mtx_test)

# -- 3. Metric plots to facilitate for additional filtering --
MetricPlots(rhapsody_obj)

# -- 4. Apply the whitelisted CellTags generated from assessment --
# Read the RDS file and get the object
rhapsody_obj <- MatrixFilteringByWhitelist(rhapsody_obj, "./rhapsody_whitelist_Mixed_99.csv")

# -- 5. Check metric plots after whitelist filtering --
MetricPlots(rhapsody_obj)

# -- 6. Additional filtering --
# Filter out cells with less than 2 CellTags ###Changed this###
rhapsody_obj <- MetricBasedFiltering(rhapsody_obj, 2, comparison = "less")

# -- 7. Last check of metric plots --
MetricPlots(rhapsody_obj)

# -- 8. Clone Calling --
# --- I. Jaccard Analysis ---
# This calculates pairwise Jaccard similarities among cells using the filtered CellTag UMI count matrix. This function takes the CellTag object with metric filtering carried out. This will generate a Jaccard similarity matrix, which is saved as a part of the object in a slot - "jaccard_mtx". It also plots a correlation heatmap with cells ordered by hierarchical clustering.
rhapsody_obj <- JaccardAnalysis(rhapsody_obj)

# --- II. Clone Calling ---
# Based on the Jaccard similarity matrix, we can call clones of cells. A clone will be selected if the correlations inside of the clones passes the cutoff given (here, 0.7 is used. It can be changed based on the heatmap/correlation matrix generated above). Using this part, a list containing the clonal identities of all cells and the count information for each clone will be stored in the object in slots - "clone.composition" and "clone.size.info".

# Clonal Identity Table clone.composition
# clone.id  | cell.barcode
# ------------------------
# Clonal ID | Cell BC

# Count Table clone.size.info
# Clone.ID  | Frequency
# ---------------------------------------
# Clonal ID | The cell number in the clone

# Call clones
rhapsody_obj <- CloneCalling(rhapsody_obj, correlation_cutoff=0.99) 

# Check them out!
rhapsody_obj@clone_composition[["v3"]]
rhapsody_obj@clone_sizes[["v3"]]

clones_df_Mixed<-rhapsody_obj@clone_composition[["v3"]]

### Save and repeat for each replicate
saveRDS(clones_df_Mixed, file = "Mixed_Clones")
saveRDS(clones_df_Rep1, file = "Mixed_Clones")
saveRDS(clones_df_Rep2, file = "Mixed_Clones")

### Add tags and combine
clones_df_Mixed$Clone_ID <- paste("Rep1_F", clones_df_Mixed$Clone_ID, sep = "_")
clones_df_Mixed$Barcode <- paste("Rep1_F", clones_df_Mixed$Barcode,  sep = "_")
clones_df_Rep1$Clone_ID <- paste("Rep2_M", clones_df_Rep1$Clone_ID, sep = "_")
clones_df_Rep1$Barcode <- paste("Rep2_M", clones_df_Rep1$Barcode,  sep = "_")
clones_df_Rep2$Clone_ID <- paste("Mixed", clones_df_Rep2$Clone_ID, sep = "_")
clones_df_Rep2$Barcode <- paste("Mixed", clones_df_Rep2$Barcode,  sep = "_")

All_Clones <- rbind(clones_df_Mixed, clones_df_Rep1, clones_df_Rep2)



###Filter Clones for clones with members of unknown clusters
#New Obj
P4_Obj_Clones <- P4_Clean
#extracting data
All_Clones_data <- All_Clones
row.names(All_Clones_data) <- All_Clones$Barcode
#Add clones
P4_Obj_Clones <- AddMetaData(object = P4_Obj_Clones, metadata = All_Clones_data)

#Plots
DimPlot(P4_Obj_Clones, reduction = "umap", label = FALSE, pt.size = 1, sizes.highlight = 1, group.by = "Clone_ID", na.value = "Grey", order = "Clone_ID">0) + NoLegend()

#Isolating only cells which are in a clone
P4_Only_Clones <- subset(P4_Obj_Clones, subset = Barcode %in% All_Clones$Barcode)
bars<-P4_Only_Clones$Barcode
res<-P4_Only_Clones$seurat_clusters
Cell_Indent_df<-data.frame(res, bars)
All_Clones_Clusters <- merge(All_Clones, Cell_Indent_df, by.x="Barcode", by.y="bars", all=TRUE)

#Identifying cells without a cluster_ID 
All_Clones_Na_Clusters <- All_Clones_Clusters[is.na(All_Clones_Clusters$res) > 0,]
NA_CLONE_IDs <- unique(All_Clones_Na_Clusters$Clone_ID)
NA_CLONE_IDs_without_triplets <- NA_CLONE_IDs

#Print the number of cells in each clone as some NA containing clones may have more than 2 cells
for (i in 1:10) {print(All_Clones_Clusters$Clone_ID[which(All_Clones_Clusters$Clone_ID == NA_CLONE_IDs[i])])}

#print NA clones with more than 3 cells
print(All_Clones_Clusters$res[which(All_Clones_Clusters$Clone_ID == "Rep1_F_11")])
print(All_Clones_Clusters$res[which(All_Clones_Clusters$Clone_ID == "Rep1_F_23")])
print(All_Clones_Clusters$res[which(All_Clones_Clusters$Clone_ID == "Rep1_F_43")])
print(All_Clones_Clusters$res[which(All_Clones_Clusters$Clone_ID == "Rep1_F_38")])
print(All_Clones_Clusters$res[which(All_Clones_Clusters$Clone_ID == "Rep1_F_115")])

#Remove clones with 2 or more remaining cell fates after removing NA clones
NA_CLONE_IDs_without_triplets <- NA_CLONE_IDs_without_triplets[NA_CLONE_IDs_without_triplets != "Rep1_F_11"]
NA_CLONE_IDs_without_triplets <- NA_CLONE_IDs_without_triplets[NA_CLONE_IDs_without_triplets != "Rep1_F_23"]
NA_CLONE_IDs_without_triplets <- NA_CLONE_IDs_without_triplets[NA_CLONE_IDs_without_triplets != "Rep1_F_43"]
NA_CLONE_IDs_without_triplets <- NA_CLONE_IDs_without_triplets[NA_CLONE_IDs_without_triplets != "Rep1_F_38"]
NA_CLONE_IDs_without_triplets <- NA_CLONE_IDs_without_triplets[NA_CLONE_IDs_without_triplets != "Rep1_F_115"]

NA_CLONE_IDs_Remove <- NA_CLONE_IDs_without_triplets

for (i in 1:5) {print(which(All_Clones_data$Clone_ID == NA_CLONE_IDs_Remove[i]))}

NA_Remove <- c(11, 12, 154, 155, 163, 164, 257, 258, 364, 365)

Keep_Clones <- setdiff(1:546, NA_Remove)

All_Good_Clones <- All_Clones[Keep_Clones]

saveRDS(All_Good_Clones, "/example_path/All_Good_Clones.RDS")

All_Good_Clones<-readRDS("All_Good_Clones_23Oct2024_99_99.RDS")

cloneIDs<-unique(All_Good_Clones$Clone_ID)

#216 unique clones, 531 cells after removing NA cells

### ##  ##     ##  #  ##  ### #### #### 
# # # # # #   #    # #  # # # ###  #
### # # # #   #    # #  # # # #     ##
# # ##  ##     ##  #  ##  # # ###  ####

#New obj
P4_Obj_Clones <- P4_Clean

#Adding clones
All_Good_Clones_data <- Clones_For_Export
row.names(All_Good_Clones_data) <- All_Good_Clones$Barcode
P4_Obj_Clones <- AddMetaData(object = P4_Obj_Clones, metadata = All_Good_Clones_data)

Clonelist_For_Export <- readRDS("All_Good_Clones_23Oct2024_99_99.RDS")
Clones_For_Export <- data.frame(1:531)
Clones_For_Export$ClusterID <- Clonelist_For_Export$seurat_clusters
Clones_For_Export$CloneID <- Clonelist_For_Export$Clone_ID
Clones_For_Export$Barcode <- Clonelist_For_Export$Barcode
Clone_Names<-unique(Clones_For_Export$CloneID)

#Check
NewCloneList

#Plot
DimPlot(P4_Obj_Clones, reduction = "umap", label = FALSE, pt.size = 1, sizes.highlight = 2, group.by = "ClusterID", na.value = "Grey", order = "Clone_ID">0) + NoLegend()

#Save
saveRDS(P4_Obj_Clones, file = "/example_path/P4_pbacBarcode_Obj.RDS")
P4_pbacBarcode_Obj <- readRDS("P4_pbacBarcode_Obj.RDS")

# # # # # # This is the Final P4 pbacBarcode object and is used from here on out for all analysis # # # # # #

#Plot
DimPlot(P4_pbacBarcode_Obj, reduction = "umap", label = FALSE, pt.size = 0.5, sizes.highlight = 1, cells.highlight = P4_Clone_Identities_Highlights, na.value = "grey")



##    ####    ####      #     #   #
# #   #  #    #        # #    #  #
#  #  # #     #       #   #   # #
###   ##      ####    #####   ##
#  #  # #     #       #   #   # #
# #   #  #    #       #   #   #  #
##    #   #   ####    #   #   #   #



##### SLINGSHOT Analysis #####
#Setting up a palette
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
dimred <- P4_pbacBarcode_Obj@reductions$umap@cell.embeddings
clustering <- P4_pbacBarcode_Obj$SCT_snn_res.0.8
counts <- as.matrix(P4_pbacBarcode_Obj@assays$RNA@counts[P4_pbacBarcode_Obj@assays$RNA@var.features, ])
lineages <- getLineages(data = dimred, clusterLabels = clustering,  start.clus = c("0", "9", "11"))
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
curves <- getCurves(lineages, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
sling_dimred<-slingshot(dimred, clustering, reducedDim = "PCA", start.clus = "0", stretch = 0, omega=TRUE)

#Plot Slingshot
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(SlingshotDataSet(sling_dimred), lwd=2, col='black', type = "curves", show.constraints = 'True')

##    ####    ####      #     #   #
# #   #  #    #        # #    #  #
#  #  # #     #       #   #   # #
###   ##      ####    #####   ##
#  #  # #     #       #   #   # #
# #   #  #    #       #   #   #  #
##    #   #   ####    #   #   #   #

##### Monocle3 analysis #####
Clean_Mono <- as.cell_data_set(P4_pbacBarcode_Obj)
Clean_Mono <- cluster_cells(Clean_Mono)
p1_Mon <- plot_cells(Clean_Mono, show_trajectory_graph = FALSE)
p2_Mon <- plot_cells(Clean_Mono, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1_Mon, p2_Mon)
Clean_Mono <- learn_graph(Clean_Mono, use_partition = TRUE)

#plot
plot_cells(Clean_Mono, color_cells_by = 'cluster', label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 4)

#Set Root Cells by selecting RGPs 
Clean_Mono <- order_cells(Clean_Mono)

#plot
plot_cells(Clean_Mono, color_cells_by = 'pseudotime', label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 4)

##### Generating a new UMAP using Monocle3's pipeline #####
IUE_Data_Mono_1<-read.csv("Mixed.csv", header=TRUE, comment.char = "#", row.names = 1)
IUE_Data_Mono_1<-as.matrix(t(IUE_Data_Mono_1))
genedata_Mono_1 <- as.data.frame(x = row.names(IUE_Data_Mono_1), row.names = row.names(IUE_Data_Mono_1))
colnames(genedata_Mono_1) <- "gene_short_name"
IUE_Data_Mono_1 <- new_cell_data_set(IUE_Data_Mono_1, gene_metadata = genedata_Mono_1)

IUE_Data_Mono_2<-read.csv("Replicate1.csv", header=TRUE, comment.char = "#", row.names = 1)
IUE_Data_Mono_2<-as.matrix(t(IUE_Data_Mono_2))
genedata_Mono_2 <- as.data.frame(x = row.names(IUE_Data_Mono_2), row.names = row.names(IUE_Data_Mono_2))
colnames(genedata_Mono_2) <- "gene_short_name"
IUE_Data_Mono_2 <- new_cell_data_set(IUE_Data_Mono_2, gene_metadata = genedata_Mono_2)

IUE_Data_Mono_3<-read.csv("Replicate2.csv", header=TRUE, comment.char = "#", row.names = 1)
IUE_Data_Mono_3<-as.matrix(t(IUE_Data_Mono_3))
genedata_Mono_3 <- as.data.frame(x = row.names(IUE_Data_Mono_3), row.names = row.names(IUE_Data_Mono_3))
colnames(genedata_Mono_3) <- "gene_short_name"
IUE_Data_Mono_3 <- new_cell_data_set(IUE_Data_Mono_3, gene_metadata = genedata_Mono_3)

cds_list <- list(IUE_Data_Mono_1, IUE_Data_Mono_2, IUE_Data_Mono_3)
New_Mono <- combine_cds(cds_list = cds_list)
New_Mono <- preprocess_cds(New_Mono, method = "PCA", num_dim = 150)
New_Mono <- reduce_dimension(New_Mono, reduction_method = "UMAP")
New_Mono <- cluster_cells(New_Mono)
p1_Mono <- plot_cells(New_Mono, show_trajectory_graph = FALSE)
p2_Mono <- plot_cells(New_Mono, color_cells_by = "partition", show_trajectory_graph = FALSE)
#plot
wrap_plots(p1_Mono, p2_Mono)
plot_cells(New_Mono, reduction_method = c("UMAP"), group_cells_by = "cluster", genes = c("Birc5", "Exo1", "Aurkb", "Gsx2", "Ascl1", "Shmt1"))
plot_cells(New_Mono, reduction_method = c("UMAP"), group_cells_by = "cluster", genes = c("Aldh1l1", "Cdh2", "Pax6", "Vim", "Nes", "Hes1", "Hes3", "Hes5", "Slc1a3", "Sox2", "Sox8"))

#Learn graph
New_Mono <- learn_graph(New_Mono, use_partition = FALSE)
#plot
plot_cells(New_Mono, color_cells_by = 'cluster', label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 4)

#Set root cells by selecting RGPs
New_Mono <- order_cells(New_Mono)

#plot
plot_cells(New_Mono, color_cells_by = 'pseudotime', label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 4)

##    ####    ####      #     #   #
# #   #  #    #        # #    #  #
#  #  # #     #       #   #   # #
###   ##      ####    #####   ##
#  #  # #     #       #   #   # #
# #   #  #    #       #   #   #  #
##    #   #   ####    #   #   #   #

##### RaceID/StemID analysis #####
RaceID_Obj  <- RaceID::Seurat2SCseq(P4_pbacBarcode_Obj)
RaceID_expData <- RaceID::getExpData(RaceID_Obj)
RaceID_res <- RaceID::pruneKnn(RaceID_expData,no_cores=10)
cl2 <- RaceID::graphCluster(RaceID_res,pvalue=0.01)
cl2$fr$V1 <- RaceID_Obj@umap[["UMAP_1"]]
cl2$fr$V2 <- RaceID_Obj@umap[["UMAP_2"]]
RaceID_Obj <- RaceID::updateSC(RaceID_Obj,res=RaceID_res,cl=cl2)
probs <- RaceID::transitionProbs(RaceID_res,cl2,pvalue=0.01)
RaceID::plotTrProbs(RaceID_Obj,probs,um=TRUE, tp =1, cex = 0.5)

##    ####    ####      #     #   #
# #   #  #    #        # #    #  #
#  #  # #     #       #   #   # #
###   ##      ####    #####   ##
#  #  # #     #       #   #   # #
# #   #  #    #       #   #   #  #
##    #   #   ####    #   #   #   #


##### P4_hdWGCNA #####
###WGCNA GO analysis via ClusterProfiler
P4_WGCNA_Clean <- P4_pbacBarcode_Obj
P4_WGCNA_Clean <- SetupForWGCNA(P4_WGCNA_Clean, gene_select = "variable", fraction = 0.05, wgcna_name = "hdWGCNA-review-test")
P4_WGCNA_Clean <- MetacellsByGroups(seurat_obj = P4_WGCNA_Clean, group.by = c("orig.ident", "seurat_clusters"), reduction = 'umap', k = 25, max_shared = 10, ident.group = 'seurat_clusters')
P4_WGCNA_Clean <- NormalizeMetacells(P4_WGCNA_Clean)

P4_WGCNA_Clean <- SetDatExpr(P4_WGCNA_Clean,group_name = "IUE_Data_Rep1", group.by='orig.ident', assay = 'SCT', slot = 'data')
P4_WGCNA_Clean <- TestSoftPowers(P4_WGCNA_Clean, networkType = 'signed')

# plot
plot_list <- PlotSoftPowers(P4_WGCNA_Clean)
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(P4_WGCNA_Clean)

# construct co-expression network:
P4_WGCNA_Clean <- ConstructNetwork(P4_WGCNA_Clean, tom_name = 'IUE_Data_Rep1')
PlotDendrogram(P4_WGCNA_Clean, main='P4 hdWGCNA Dendrogram')
TOM <- GetTOM(P4_WGCNA_Clean)

# generating modules
P4_WGCNA_Clean <- ModuleEigengenes(P4_WGCNA_Clean,group.by.vars="orig.ident")
hMEs <- GetMEs(P4_WGCNA_Clean)
MEs <- GetMEs(P4_WGCNA_Clean, harmonized=FALSE)
P4_WGCNA_Clean <- ModuleConnectivity(P4_WGCNA_Clean, group.by = 'orig.ident', group_name = 'IUE_Data_Rep1')
P4_WGCNA_Clean <- ResetModuleNames(P4_WGCNA_Clean, new_name = "Rep1_")

# plot
p <- PlotKMEs(P4_WGCNA_Clean, ncol=5)

# module assignment table:
modules <- GetModules(P4_WGCNA_Clean) %>% subset(module != 'grey')
hub_df <- GetHubGenes(P4_WGCNA_Clean, n_hubs = 10)

library(UCell)
P4_WGCNA_Clean <- ModuleExprScore(P4_WGCNA_Clean, n_genes = 25, method='UCell')

# plot
plot_list <- ModuleFeaturePlot(P4_WGCNA_Clean, features='hMEs', order=TRUE)
wrap_plots(plot_list, ncol=6)

plot_list <- ModuleFeaturePlot(P4_WGCNA_Clean, features='scores', order='shuffle', ucell = TRUE)
wrap_plots(plot_list, ncol=6)

# Final modules
hdTurquoisemodule <- modules[order(modules$kME_Rep1_1 , decreasing = TRUE),]
hdTurquoiseGenes <- unlist2(hdTurquoisemodule$gene_name[which(hdTurquoisemodule$module == "Rep1_1")])

hdBrownmodule <- modules[order(modules$kME_Rep1_2 , decreasing = TRUE),]
hdBrownGenes <- unlist2(hdBrownmodule$gene_name[which(hdBrownmodule$module == "Rep1_2")])

hdBluemodule <- modules[order(modules$kME_Rep1_3 , decreasing = TRUE),]
hdBlueGenes <- unlist2(hdBluemodule$gene_name[which(hdBluemodule$module == "Rep1_3")])

hdGreenmodule <- modules[order(modules$kME_Rep1_4 , decreasing = TRUE),]
hdGreenGenes <- unlist2(hdGreenmodule$gene_name[which(hdGreenmodule$module == "Rep1_4")])

hdYellowmodule <- modules[order(modules$kME_Rep1_5 , decreasing = TRUE),]
hdYellowGenes <- unlist2(hdYellowmodule$gene_name[which(hdYellowmodule$module == "Rep1_5")])

# plot
FeaturePlot(P4_Clean, features = c(hdTurquoiseGenes[1:12]))
FeaturePlot(P4_Clean, features = c(hdBrownGenes[1:12]))
FeaturePlot(P4_Clean, features = c(hdBlueGenes[1:12]))
FeaturePlot(P4_Clean, features = c(hdGreenGenes[1:12]))
FeaturePlot(P4_Clean, features = c(hdYellowGenes[1:12]))

#GO analysis
organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
keytypes(org.Mm.eg.db)
yellow.GO.enriched <- enrichGO(hdYellowGenes, org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)
require(DOSE)
dotplot(yellow.GO.enriched, showCategory=20)



##    ####    ####      #     #   #
# #   #  #    #        # #    #  #
#  #  # #     #       #   #   # #
###   ##      ####    #####   ##
#  #  # #     #       #   #   # #
# #   #  #    #       #   #   #  #
##    #   #   ####    #   #   #   #



##############################################################
#    ###   ####      #   #########  #####     #    #        #
#   #      #   #    # #      #        #      # #   #        #
#   #      #   #   #   #     #        #     #   #  #        #
#    #     #   #   #   #     #        #     #   #  #        #
#     #    ####    #####     #        #     #####  #        #
#      #   #       #   #     #        #     #   #  #        #
#      #   #       #   #     #        #     #   #  #        #
#     #    #       #   #     #        #     #   #  #        #
#  ###     #       #   #     #      #####   #   #  #######  #
##############################################################

# The original P4 and P22 Curio objects were generated by Curio's analysis pipeline on DNAnexus and were imported into R
# The output of the Curio pipeline for the P4 scSpatialseq was P4_WT_Cortex_and_Ventricle_Deeper_seurat.RDS
# The output of the Curio pipeline for the P22 scSpatialseq was P22_WT_Cortex_and_Ventricle_Deeper_seurat.RDS
# We take these seurat objects and used them to generate RCTD using the following method

#P4
P4_WT_Cortex_and_Ventricle_Deeper_seurat <- readRDS("P4_WT_Cortex_and_Ventricle_Deeper_seurat.RDS")
str(P4_WT_Cortex_and_Ventricle_Deeper_seurat) #this is the P4 obj

#P22
P22_WT_Cortex_and_Ventricle_Deeper_seurat <- readRDS("P22_WT_Cortex_and_Ventricle_Deeper_seurat.RDS")
str(P22_WT_Cortex_and_Ventricle_Deeper_seurat) #this is the P22 obj

# set up reference for RCTD, we used the original P4_Obj here before we added clones and removed the unknown populations
P4_Obj_RCTD <- P4_Obj

ref <- P4_Obj_RCTD
ref <- UpdateSeuratObject(ref)
Idents(ref) <- "celltype"
str(ref)
# Generating REF from P4_Obj
Ref_counts_RCTD <- ref[["RNA"]]@counts
Ref_cluster_RCTD <- as.factor(ref$seurat_clusters)
names(Ref_cluster_RCTD) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(Ref_counts_RCTD, Ref_cluster_RCTD, nUMI)

# Generating Query from P4_WT_Cortex_and_Ventricle_Deeper_seurat
Curio_P4_RCTD_Query <- P4_WT_Cortex_and_Ventricle_Deeper_seurat
Query_counts_RCTD <- Curio_P4_RCTD_Query[["RNA"]]@counts
coords <- GetTissueCoordinates(Curio_P4_RCTD_Query)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, Query_counts_RCTD, colSums(Query_counts_RCTD))

#Run RCTD
myRCTD <- create.RCTD(query, reference, max_cores = 11, CELL_MIN_INSTANCE = 20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
Curio_P4_RCTD_Query <- AddMetaData(Curio_P4_RCTD_Query, metadata = myRCTD@results$results_df)

#save it
saveRDS(Curio_P4_RCTD_Query, file = "P4_Curio_Spatial_Obj.RDS")
P4_Curio_Spatial_Obj <- readRDS("/example_path/P4_Curio_Spatial_Obj.RDS")
saveRDS(myRCTD, file = "P4_myRCTD.RDS")
P4_myRCTD <- readRDS("/example_path/P4_myRCTD.RDS")

#plot predictions
p1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, group.by = "first_type")
p2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, group.by = "second_type")
p1 | p2

#These will highlight each cluster in a spatial plot using the cluster identities from the P4_Obj
first_type <- data.frame(P4_Curio_Spatial_Obj$first_type)
first_type_0 <- row.names(first_type)[which(first_type == 0)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_0, cols.highlight = c("red", "white"))
second_type <- data.frame(P4_Curio_Spatial_Obj$second_type)
second_type_0 <- row.names(second_type)[which(second_type == 0)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_0, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_1 <- row.names(first_type)[which(first_type == 1)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_1, cols.highlight = c("red", "white"))
second_type_1 <- row.names(second_type)[which(second_type == 1)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_1, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_2 <- row.names(first_type)[which(first_type == 2)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_2, cols.highlight = c("red", "white"))
second_type_2 <- row.names(second_type)[which(second_type == 2)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_2, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_3 <- row.names(first_type)[which(first_type == 3)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_3, cols.highlight = c("red", "white"))
second_type_3 <- row.names(second_type)[which(second_type == 3)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_3, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_4 <- row.names(first_type)[which(first_type == 4)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_4, cols.highlight = c("red", "white"))
second_type_4 <- row.names(second_type)[which(second_type == 4)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_4, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_5 <- row.names(first_type)[which(first_type == 5)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_5, cols.highlight = c("red", "white"))
second_type_5 <- row.names(second_type)[which(second_type == 5)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_5, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_6 <- row.names(first_type)[which(first_type == 6)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_6, cols.highlight = c("red", "white"))
second_type_6 <- row.names(second_type)[which(second_type == 6)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_6, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_7 <- row.names(first_type)[which(first_type == 7)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_7, cols.highlight = c("red", "white"))
second_type_7 <- row.names(second_type)[which(second_type == 7)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_7, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_8 <- row.names(first_type)[which(first_type == 8)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_8, cols.highlight = c("red", "white"))
second_type_8 <- row.names(second_type)[which(second_type == 8)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_8, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_9 <- row.names(first_type)[which(first_type == 9)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_9, cols.highlight = c("red", "white"))
second_type_9 <- row.names(second_type)[which(second_type == 9)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_9, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_10 <- row.names(first_type)[which(first_type == 10)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_10, cols.highlight = c("red", "white"))
second_type_10 <- row.names(second_type)[which(second_type == 10)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_10, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_11 <- row.names(first_type)[which(first_type == 11)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_11, cols.highlight = c("red", "white"))
second_type_11 <- row.names(second_type)[which(second_type == 11)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_11, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_12 <- row.names(first_type)[which(first_type == 12)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_12, cols.highlight = c("red", "white"))
second_type_12 <- row.names(second_type)[which(second_type == 12)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_12, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_13 <- row.names(first_type)[which(first_type == 13)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_13, cols.highlight = c("red", "white"))
second_type_13 <- row.names(second_type)[which(second_type == 13)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_13, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_14 <- row.names(first_type)[which(first_type == 14)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_14, cols.highlight = c("red", "white"))
second_type_14 <- row.names(second_type)[which(second_type == 14)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_14, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_15 <- row.names(first_type)[which(first_type == 15)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_15, cols.highlight = c("red", "white"))
second_type_15 <- row.names(second_type)[which(second_type == 15)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_15, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_16 <- row.names(first_type)[which(first_type == 16)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_16, cols.highlight = c("red", "white"))
second_type_16 <- row.names(second_type)[which(second_type == 16)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_16, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_17 <- row.names(first_type)[which(first_type == 17)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_17, cols.highlight = c("red", "white"))
second_type_17 <- row.names(second_type)[which(second_type == 17)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_17, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_18 <- row.names(first_type)[which(first_type == 18)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_18, cols.highlight = c("red", "white"))
second_type_18 <- row.names(second_type)[which(second_type == 18)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_18, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_19 <- row.names(first_type)[which(first_type == 19)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_19, cols.highlight = c("red", "white"))
second_type_19 <- row.names(second_type)[which(second_type == 19)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_19, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

first_type_20 <- row.names(first_type)[which(first_type == 20)]
Spatial_P1 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = first_type_20, cols.highlight = c("red", "white"))
second_type_20 <- row.names(second_type)[which(second_type == 20)]
Spatial_P2 <- SpatialDimPlot(P4_Curio_Spatial_Obj, cells.highlight = second_type_20, cols.highlight = c("red", "white"))
Spatial_P1 | Spatial_P2

##We repeat the process for P22 Spatial Data
# set up query
Curio_P22_RCTD_Query <- P22_WT_Cortex_and_Ventricle_Deeper_seurat
Query_counts_RCTD <- Curio_P22_RCTD_Query[["RNA"]]@counts
coords <- GetTissueCoordinates(Curio_P22_RCTD_Query)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query22 <- SpatialRNA(coords, Query_counts_RCTD, colSums(Query_counts_RCTD))

#Run RCTD
myRCTD <- create.RCTD(query22, reference, max_cores = 11, CELL_MIN_INSTANCE = 20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
Curio_P22_RCTD_Query <- AddMetaData(Curio_P22_RCTD_Query, metadata = myRCTD@results$results_df)

#Save them
saveRDS(Curio_P22_RCTD_Query, file = "P22_Curio_Spatial_Obj.RDS")
P22_Curio_Spatial_Obj <- readRDS("/example_path/P22_Curio_Spatial_Obj.RDS")
saveRDS(myRCTD, file = "P22_myRCTD.RDS")
P22_myRCTD <- readRDS("/example_path/P22_myRCTD.RDS")

#plot
p1 <- SpatialDimPlot(Curio_P22_RCTD, group.by = "first_type")
p2 <- SpatialDimPlot(Curio_P22_RCTD, group.by = "second_type")
p1 | p2

#The same first_type identity plots can be generated as with P4 by taking the previous section of code and repeating it with P22_Curio_Spatial_Obj


##    ####    ####      #     #   #
# #   #  #    #        # #    #  #
#  #  # #     #       #   #   # #
###   ##      ####    #####   ##
#  #  # #     #       #   #   # #
# #   #  #    #       #   #   #  #
##    #   #   ####    #   #   #   #



##### Printing the unique colours for our plots ####
Curio_P4_RCTD_colour <- Curio_P4_RCTD
Curio_P22_RCTD_colour <- Curio_P22_RCTD

coloursforplot2 <- print("olivedrab3", "deeppink3", "chocolate", "darkslategrey", "orchid1", "darkturquoise", "indianred2", "yellow2", "maroon", "goldenrod2", "greenyellow", "firebrick2", "cyan1", "forestgreen", "khaki2", "deeppink4", "palevioletred4", "thistle2", "dodgerblue4", "lightsteelblue1", "tan3")
my_cols <-  c("violetred4", "deeppink3", "deepskyblue", "seagreen", "indianred2", "turquoise", "purple3", "yellowgreen", "orangered1", "sienna2", "lawngreen", "navy", "goldenrod", "darkmagenta", "thistle1", "orchid4", "yellow3")

Curio_P4_RCTD_colour$first_type <- factor(Curio_P4_RCTD_colour$first_type, levels = names(my_cols))
SpatialDimPlot(Curio_P4_RCTD_colour, group.by = "first_type", cols = my_cols)

Curio_P22_RCTD_colour$first_type <- factor(Curio_P22_RCTD_colour$first_type, levels = names(my_cols))
SpatialDimPlot(Curio_P22_RCTD_colour, group.by = "first_type", cols = my_cols)

#You can do the same for second predictions
Curio_P4_RCTD_colour$second_type <- factor(Curio_P4_RCTD_colour$second_type, levels = names(my_cols))
SpatialDimPlot(Curio_P4_RCTD_colour, group.by = "second_type", cols = my_cols)

Curio_P22_RCTD_colour$second_type <- factor(Curio_P22_RCTD_colour$second_type, levels = names(my_cols))
SpatialDimPlot(Curio_P22_RCTD_colour, group.by = "second_type", cols = my_cols)

P4_pbacBarcode_Obj$colours <- P4_pbacBarcode_Obj$seurat_clusters
P4_pbacBarcode_Obj$colours <- factor(P4_pbacBarcode_Obj$colours, levels = names(my_cols))
DimPlot(P4_pbacBarcode_Obj, group.by = "colours", cols = my_cols)
