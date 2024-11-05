library(Seurat) # 
library(Matrix)
library(ggplot2)
library(sctransform)
library(stringr)
library(cowplot) # used for CCA
library(patchwork) # used for CCA
library(dplyr) # used for print markers
library(here) # notes path for easy access in commands
library(future)
options(future.globals.maxSize = 8000 * 1024^2)

#UMAPcolors <- c('palevioletred', 'palevioletred1', 'palevioletred2', 'palevioletred3',
		                'pink1', 'pink2', 'pink3',
				                'hotpink2','hotpink3',
				                'skyblue1', 'skyblue2',
						                'steelblue', 'steelblue1', 'steelblue2', 'steelblue3', 'steelblue4', 'dodgerblue3',
						                'darkseagreen1', 'darkseagreen2', 'darkseagreen',
								                'lightgoldenrod', 'lightgoldenrod1', 'lightgoldenrod2', 'gold', 'gold2', 'goldenrod1', 'goldenrod3',
								                'sienna1', 'sienna2', 'sienna3', 'tomato3',
										                'red',
										                'purple4')
#
#UMAPcolors_all <- c('gray0', 'darkgray', 'darkgray', 'darkgray', 'darkgray', 'darkgray', 'darkgray', 'darkgray',
		                        'palevioletred', 'palevioletred1', 'palevioletred2', 'palevioletred3',
					                    'pink1', 'pink2', 'pink3',
					                    'hotpink2','hotpink3',
							                        'skyblue1', 'skyblue2',
							                        'steelblue', 'steelblue1', 'steelblue2', 'steelblue3', 'steelblue4', 'dodgerblue3',
										                    '    darkseagreen1', 'darkseagreen2', 'darkseagreen',
										                    'lightgoldenrod', 'lightgoldenrod1', 'lightgoldenrod2', 'gold', 'gold2', 'goldenrod1', 'goldenrod3',
												                        'sienna1', 'sienna2', 'sienna3', 'tomato3',
												                        'red',
															                    'purple4')
#
## load whole optic lobe data
#allnorm.integrated.tree.ordered <- readRDS(here("data","Final-ALL_UMAP-ordered_1.rds"))
#
## Update Seurat object to v5, add metadata
#whole_ol<- UpdateSeuratObject(object = allnorm.integrated.tree.ordered)
#whole_ol<- AddMetaData(whole_ol, "whole_ol", col.name = "sample")
#whole_ol@project.name <- "whole_ol"
#
##### Read in topo sample CellRanger output cells and genes
#b <- Read10X_h5(here("data","topo_sample_cellranger_filtered_counts", "h5","B_filtered_feature_bc_matrix.h5"))
#t <- Read10X_h5(here("data","topo_sample_cellranger_filtered_counts", "h5","T_filtered_feature_bc_matrix.h5"))
#r <- Read10X_h5(here("data","topo_sample_cellranger_filtered_counts", "h5","R_filtered_feature_bc_matrix.h5"))
#l <- Read10X_h5(here("data","topo_sample_cellranger_filtered_counts", "h5","L_filtered_feature_bc_matrix.h5"))
#
#####Creating Seurat objects for each sample type
#bsc <- CreateSeuratObject(counts = b, project = "bottom", min.cells = 3, min.features = 500)
#bsc<- AddMetaData(bsc, "bottom", col.name = "sample")
##Idents(bsc)<- "bottom"
#tsc <- CreateSeuratObject(counts = t, project = "top", min.cells = 3, min.features = 500)
#tsc<- AddMetaData(tsc, "top", col.name = "sample")
##Idents(tsc)<- "top"
#rsc <- CreateSeuratObject(counts = r, project = "right", min.cells = 3, min.features = 500)
#rsc<- AddMetaData(rsc, "right", col.name = "sample")
##Idents(rsc)<- "right"
#lsc <- CreateSeuratObject(counts = l, project = "left", min.cells = 3, min.features = 500)
#lsc<- AddMetaData(lsc, "left", col.name = "sample")
##Idents(lsc)<- "left"
#
#topos<- merge(bsc, y = c(tsc,rsc,lsc), add.cell.ids = c("b","t","r","l"), project = "all_topos")
#
#topos
#
## Remove cells with high mitochondrial content
#Thres_nCount_topos <- VlnPlot(topos, features = "nCount_RNA", y.max = 35000, pt.size = 0) + geom_hline(yintercept=1000, linetype='dotted', col = 'black') + geom_hline(yintercept=20000, linetype='dotted', col = 'black')
#
#Thres_nFeat_topos <- VlnPlot(topos, features = "nFeature_RNA", y.max = 8000, pt.size = 0) + geom_hline(yintercept=600, linetype='dotted', col = 'black')
#
#topos[["percent.mt"]]<- PercentageFeatureSet(topos, pattern = "^[^obimac]")
#grep( "^[^obimac]", rownames(topos), value = T)
#
#Thres_MT_topos <- VlnPlot(topos, features = "percent.mt", y.max = 30, pt.size = 0) + geom_hline(yintercept=6, linetype='dotted', col = 'black') 
#
#pdf("topo_nCount_nFeat_MT.pdf")
#Thres_nCount_topos
#Thres_nFeat_topos
#Thres_MT_topos
#dev.off()
#
#topos 
#
#topos<- subset(topos, subset = percent.mt < 6 & nCount_RNA > 1000 & nCount_RNA < 20000 & 
	       	                        #nFeature_RNA > 600)
#topos
#
#topos<- SCTransform(topos)
#saveRDS(topos, here("data","topos.RDS"))

topos <- readRDS(here("data","topos.RDS"))
  
topos <- RunPCA(topos, verbose = FALSE)
topos <- RunUMAP(topos, dims = 1:30, verbose = FALSE)
topos <- RunTSNE(topos, dims = 1:30, verbose = FALSE)

saveRDS(topos, here("data","topos.RDS"))

topos <- FindNeighbors(topos, dims = 1:30, verbose = FALSE)
topos <- FindClusters(topos, verbose = FALSE)
DimPlot(topos, label = TRUE)



topos_deg <- FindAllMarkers(topos, only.pos = TRUE)

write_tsv(topos_deg, here("data","20241010_topos_deg.tsv"))






#
##anchors<- FindTransferAnchors(reference = whole_ol, query = topos, normalization.method = "SCT",
##			          features = VariableFeatures(object = whole_ol), reference.assay = "RNA", query.assay = "RNA", 
##				      reduction = "pcaproject")
##saveRDS(anchors, here("data","whole_ol_anchors.RDS"))
#
##predictions<- TransferData(anchorset = anchors, refdata = whole_ol$seurat_clusters)
##topos <- AddMetaData(object = topos, metadata = predictions)
#
