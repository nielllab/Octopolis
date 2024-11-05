
library(Seurat)
library(patchwork)
library(tidyverse)
library(here) #helps create platform independent file paths

#### Read in CellRanger output cells and genes

b <- Read10X(here("data","sc_OL_topo_analysis", "sc_topo_B_outs", "filtered_feature_bc_matrix"))
t <- Read10X(here("data","sc_OL_topo_analysis", "sc_topo_T_outs", "filtered_feature_bc_matrix"))
r <- Read10X(here("data","sc_OL_topo_analysis", "sc_topo_R_outs", "filtered_feature_bc_matrix"))
l <- Read10X(here("data","sc_OL_topo_analysis", "sc_topo_L_outs", "filtered_feature_bc_matrix"))


####Creating Seurat objects for each sample type

bsc <- CreateSeuratObject(counts = b, project = "bottom", min.cells = 3, min.features = 200)
bsc<- AddMetaData(bsc, "bottom", col.name = "sample")

tsc <- CreateSeuratObject(counts = t, project = "top", min.cells = 3, min.features = 200)
tsc<- AddMetaData(tsc, "top", col.name = "sample")

rsc <- CreateSeuratObject(counts = r, project = "right", min.cells = 3, min.features = 200)
rsc<- AddMetaData(rsc, "right", col.name = "sample")

lsc <- CreateSeuratObject(counts = l, project = "left", min.cells = 3, min.features = 200)
lsc<- AddMetaData(lsc, "left", col.name = "sample")


ol_sc<- readRDS(here("data","Final-ALL_UMAP-ordered_1.rds"))
whole_ol_sc<- UpdateSeuratObject(object = ol_sc)
#Layers(whole_ol_sc[["RNA"]])
whole_ol_sc<- AddMetaData(whole_ol_sc, "whole_ol", col.name = "sample")
whole_ol_sc<- 
  DefaultAssay(whole_ol_sc) <- "RNA"
whole_ol_sc[['integrated']] <- NULL
whole_ol_sc[['SCT']] <- NULL
whole_ol_sc


# create a Seurat object initialized with multiple layers
#ol_sc_topo<- merge(whole_ol_sc, y= c(bsc,tsc,rsc,lsc), add.cell.ids = c("whole_ol","bottom","top","right","left"), project = "topo")
ol_sc_topo<- merge(bsc, y= c(tsc,rsc,lsc), add.cell.ids = c("bottom","top","right","left"), project = "topo")
Layers(ol_sc_topo[["RNA"]])

ol_sc_topo

#### run standard anlaysis workflow on non-integrated data

ol_sc_topo <- NormalizeData(ol_sc_topo)
ol_sc_topo <- FindVariableFeatures(ol_sc_topo)
ol_sc_topo <- ScaleData(ol_sc_topo)
ol_sc_topo <- RunPCA(ol_sc_topo)
ol_sc_topo <- FindNeighbors(ol_sc_topo, dims = 1:30, reduction = "pca")
ol_sc_topo <- FindClusters(ol_sc_topo, resolution = 1, cluster.name = "unintegrated_clusters")

ol_sc_topo <- RunUMAP(ol_sc_topo, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
unintegrated_dimplot<- DimPlot(ol_sc_topo, reduction = "umap.unintegrated", group.by = c("sample", "seurat_clusters"))


#### Run integration on layers
ol_sc_topo <- IntegrateLayers(object = ol_sc_topo, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
ol_sc_topo[["RNA"]] <- JoinLayers(ol_sc_topo[["RNA"]])

ol_sc_topo <- FindNeighbors(ol_sc_topo, reduction = "integrated.cca", dims = 1:30)
ol_sc_topo <- FindClusters(ol_sc_topo, resolution = 1)

ol_sc_topo <- RunUMAP(ol_sc_topo, dims = 1:30, reduction = "integrated.cca")

# Visualization
#DimPlot(ol_sc_topo, reduction = "umap", group.by = c("stim", "seurat_annotations"))

#DimPlot(ol_sc_topo, reduction = "umap", split.by = "stim")

automagic::make_deps_file()
sessionInfo()
