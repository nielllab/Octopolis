# Process data and normalize before merging to see if UMAP is able to integrate all data appropriately

# Set working directory for entire run and load all libraries

## To run entire code as single command, use the Source tab

setwd("/Users/Hannah Bishop/Desktop/MeaS/Seurat")

library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(sctransform)
library(stringr)
library(dplyr)

# Data Set 1: Project "OctoSeq1_MTnorm3"
raw1 <- Read10X(data.dir = "Y:/OctoSeq_For_Seurat/OctoSeq1/raw_feature_bc_matrix")

ref <- read.csv("C:/Users/Hannah\ Bishop/Desktop/MeaS/miniNCBI_annotations_alltomitochondria_022420.csv", stringsAsFactors=FALSE)

ngenes <- length(raw1@Dimnames[[1]])
for (g in 1:ngenes){
  gene<-raw1@Dimnames[[1]][g]
  gene<-substr(gene,6,str_length(gene)-2)
  ind<-grep(gene,ref[[1]])
  if (length(ind)>0) {
    id <- ref[[ind[1],2]]
    if (str_length(id)>0) {
      id <- str_remove_all(id,"\\(") # parentheses mess up gene names as dimensions
      id <- str_remove_all(id,"\\)")
      id <- substr(id,1,60) # keep it short
      raw1@Dimnames[[1]][g]<- paste(id,gene,sep='-')
    }
  }
}

raw1 <- CreateSeuratObject(counts = raw1, project = "OctoSeq1_MTnorm3", min.cells = 3, min.features = 200)
raw1

mito.genes <- grep(pattern = "^mitochondria-", x = rownames(x = raw1), value = TRUE)
percent.mito <- Matrix::colSums(raw1) / Matrix::colSums(raw1)
raw1[["percent.mt"]] <- PercentageFeatureSet(raw1, pattern = "^mitochondria-")
raw1 <- subset(raw1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
plot1 <- FeatureScatter(raw1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(raw1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

raw1 <- NormalizeData(raw1, normalization.method = "LogNormalize", scale.factor = 10000)

raw1 <- FindVariableFeatures(raw1, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(raw1), 10)

plot3 <- VariableFeaturePlot(raw1)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))

#raw1.genes <- rownames(raw1)
#raw1.norm <- ScaleData(raw1, vars.to.regress = c("nCount_RNA", "percent.mt"))

saveRDS(raw1, file = "/Users/Hannah\ Bishop/Desktop/MeaS/Seurat/OSr1MTnorm3.rds")

# Will need to run through RunPCA command




# Data Set 2: Project "OctoSeq2_MTnorm3"

raw2 <- Read10X(data.dir = "Y:/OctoSeq_For_Seurat/OctoSeq2.1/raw_feature_bc_matrix")

ref <- read.csv("C:/Users/Hannah\ Bishop/Desktop/MeaS/miniNCBI_annotations_alltomitochondria_022420.csv", stringsAsFactors=FALSE)

ngenes <- length(raw2@Dimnames[[1]])
for (g in 1:ngenes){
  gene<-raw2@Dimnames[[1]][g]
  gene<-substr(gene,6,str_length(gene)-2)
  ind<-grep(gene,ref[[1]])
  if (length(ind)>0) {
    id <- ref[[ind[1],2]]
    if (str_length(id)>0) {
      id <- str_remove_all(id,"\\(") # parentheses mess up gene names as dimensions
      id <- str_remove_all(id,"\\)")
      id <- substr(id,1,60) # keep it short
      raw2@Dimnames[[1]][g]<- paste(id,gene,sep='-')
    }
  }
}

raw2 <- CreateSeuratObject(counts = raw2, project = "OctoSeq2_MTnorm3", min.cells = 3, min.features = 200)
raw2

mito.genes <- grep(pattern = "^mitochondria-", x = rownames(x = raw2), value = TRUE)
percent.mito <- Matrix::colSums(raw2) / Matrix::colSums(raw2)
raw2[["percent.mt"]] <- PercentageFeatureSet(raw2, pattern = "^mitochondria-")
raw2 <- subset(raw2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
plot1 <- FeatureScatter(raw2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(raw2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

raw2 <- NormalizeData(raw2, normalization.method = "LogNormalize", scale.factor = 10000)

raw2 <- FindVariableFeatures(raw2, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(raw2), 10)

plot3 <- VariableFeaturePlot(raw2)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))

#raw2.genes <- rownames(raw2)
#raw2.norm <- ScaleData(raw2, vars.to.regress = c("nCount_RNA", "percent.mt"))

saveRDS(raw2, file = "/Users/Hannah\ Bishop/Desktop/MeaS/Seurat/OSr2MTnorm3.rds")

# Will need to run through RunPCA command





# Data Set 3: Project "OctoSeq3_MTnorm3"

raw3 <- Read10X(data.dir = "Y:/OctoSeq_For_Seurat/OctoSeq2.2/raw_feature_bc_matrix")

ref <- read.csv("C:/Users/Hannah\ Bishop/Desktop/MeaS/miniNCBI_annotations_alltomitochondria_022420.csv", stringsAsFactors=FALSE)

ngenes <- length(raw3@Dimnames[[1]])
for (g in 1:ngenes){
  gene<-raw3@Dimnames[[1]][g]
  gene<-substr(gene,6,str_length(gene)-2)
  ind<-grep(gene,ref[[1]])
  if (length(ind)>0) {
    id <- ref[[ind[1],2]]
    if (str_length(id)>0) {
      id <- str_remove_all(id,"\\(") # parentheses mess up gene names as dimensions
      id <- str_remove_all(id,"\\)")
      id <- substr(id,1,60) # keep it short
      raw3@Dimnames[[1]][g]<- paste(id,gene,sep='-')
    }
  }
}

raw3 <- CreateSeuratObject(counts = raw3, project = "OctoSeq3_MTnorm3", min.cells = 3, min.features = 200)
raw3

mito.genes <- grep(pattern = "^mitochondria-", x = rownames(x = raw3), value = TRUE)
percent.mito <- Matrix::colSums(raw3) / Matrix::colSums(raw3)
raw3[["percent.mt"]] <- PercentageFeatureSet(raw3, pattern = "^mitochondria-")
raw3 <- subset(raw3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
plot1 <- FeatureScatter(raw3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(raw3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

raw3 <- NormalizeData(raw3, normalization.method = "LogNormalize", scale.factor = 10000)

raw3 <- FindVariableFeatures(raw3, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(raw3), 10)

plot3 <- VariableFeaturePlot(raw3)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))

#raw3.genes <- rownames(raw3)
#raw3.norm <- ScaleData(raw3, vars.to.regress = c("nCount_RNA", "percent.mt"))

saveRDS(raw3, file = "/Users/Hannah\ Bishop/Desktop/MeaS/Seurat/OSr3MTnorm3.rds")

# Will need to run through RunPCA command

## files saved out as all.norm - taken through pre-processing steps through to ScaleDatafunction


# Begin PCA steps, save out just in case (scale data prior)
#OSr1MT <- ScaleData(raw1.norm, vars.to.regress = c("nCount_RNA", "percent.mt"))
#OSr2MT <- ScaleData(raw2.norm, vars.to.regress = c("nCount_RNA", "percent.mt"))
#OSr3MT <- ScaleData(raw3.norm, vars.to.regress = c("nCount_RNA", "percent.mt"))

allOSnorm <- merge(x = raw1, y = c(raw2, raw3), add.cell.ids = c("run1", "run2", "run3"), project = "OSallNormMT2", merge.data = TRUE)
GetAssayData(allOSnorm)[1:10, 1:15]

allOSnorm <- FindVariableFeatures(allOSnorm)
allOSnorm <- ScaleData(allOSnorm, vars.to.regress = c("nCount_RNA", "percent.mt"))
all.norm <- RunPCA(allOSnorm, features = VariableFeatures(object = allOSnorm), npcs = 200)
#all.norm <- saveRDS(all.norm, file = "/Users/Hannah\ Bishop/Desktop/MeaS/Seurat/allOSnorm_merged2.rds")

# Follow appropriate steps to prepare for UMAP visualization

all.norm <- FindNeighbors(all.norm, dims = 1:50)
all.norm <- FindClusters(all.norm, reduction.type = "pca", dims = 1:50, resolution = 1)
all.norm <- RunUMAP(all.norm, dims = 1:50) 
DimPlot(all.norm, reduction = "umap", label = TRUE)

# Visualize Runs superimposed on UMAP

DimPlot(all.norm, reduction = "umap", group.by = "orig.ident")
#DimPlot(all.norm, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) 

# save files

all.norm <- saveRDS(all.norm, file = "/Users/Hannah\ Bishop/Desktop/MeaS/Seurat/allOSnorm_merged3.rds")

