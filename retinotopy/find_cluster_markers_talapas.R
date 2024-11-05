---
  title: "FindMarkers for each cluster vs all"
#output: html_notebook
#output: html_notebook: default
#html_document: default
pdf_document: default
---
  
  ### Load libraries

```{r, require = T}
library(Seurat) # 
library(Matrix)
library(ggplot2)
library(sctransform)
library(stringr)
library(cowplot) # used for CCA
library(patchwork) # used for CCA
library(dplyr) # used for print markers
library(here) # notes path for easy access in commands
library(scCustomize)
options("SCpubr.verbose" = FALSE)
library(SCpubr)
library(tidyverse)
library(ggplot2)
library(future)
options(future.globals.maxSize = 8000 * 1024^2)
plan("multisession")
plan()
```



fil_sct_clust_findmarkers <- readRDS(here("data","FindMarkers_filtered_SCTransform_clustered_topos.RDS"))
fil_sct_clust_findmarkers

gene_master <- read_csv(here("data","GeneIDs.csv"))



# function to add gene names to data table
# assumes first column of table is G gene
# uses master list with gene name in 1st column, G gene in 3rd column
addGeneNamesTable<- function(data,geneID, gene, desc){ # a little redundant to send in geneID but this solves the problem that different objects name it differently
  data$GeneName <- '' # create empty column
  for (i in 1:nrow(data)){
    #ggene = data$gene[[i]] #cluster markers
    obgene = geneID[[i]] # individual cluster markers
    #ggene = data$Feature[[i]] #dimheatmap
    #ggene = data$feature[[i]] #viz dim loadings
    loc <- which(gene ==obgene)
    if (length(loc)>0) {  # if something is found
      data$GeneName[[i]] <- desc[loc[1]] # if multiple found, use 1st one
    }
  }
  
  data
}

replaceLabelsX <- function(p,gene,desc){
  pp<- ggplot_build(p)
  oldxlabels = ggplot_build(p)$layout$panel_params[[1]]$x$breaks
  newxlabels = oldxlabels
  for (i in 1:length(oldxlabels)){
    loc<-which(gene==(oldxlabels[[i]]))
    if (length(loc)>0){
      newxlabels[[i]] = paste(oldxlabels[[i]], "-", desc[[loc[1]]])
      newxlabels[[i]] = substr(newxlabels[[i]],1,100)}}
  p<- p + scale_x_discrete(labels = newxlabels)
}

replaceLabelsY <- function(p,gene,desc){
  pp<- ggplot_build(p)
  oldylabels = ggplot_build(p)$layout$panel_params[[1]]$y$breaks
  newylabels = oldylabels
  for (i in 1:length(oldylabels)){
    loc<-which(gene==(oldylabels[[i]]))
    if (length(loc)> 0) {
      newylabels[[i]] = paste(oldylabels[[i]], "-", desc[[loc[1]]])
      newylabels[[i]] = substr(newylabels[[i]],1,100)}}
  p<- p + scale_y_discrete(labels = newylabels)
}

shorten <- function(x)
{
  stringr::str_trunc(x, 100)
}

# Take all cells in cluster 2, and find markers that separate cells in the 'g1' group (metadata
# variable 'group')
#markers <- FindMarkers(fil_sct_clust_findmarkers, ident.1 = "bottom", group.by = 'sample', subset.ident = "0")

# make list of cluster identities
clusters <- levels(fil_sct_clust@active.ident)

# create dotplots of top 30 markers for each cluster
png(here("output","cluster_DEG.png"), width = 12, height = 7, units = "in", res = 300)

for (m in 0:3){ #length(allnorm.integrated.tree.ordered)
  clustermarkers_i <- FindMarkers(fil_sct_clust_findmarkers, ident.1 = m, logfc.threshold = 0.25, only.pos = TRUE) #, test.use = "roc"
  clustermarkers_i$diff_exp <- (clustermarkers_i$pct.1/clustermarkers_i$pct.2)
  clustermarkers_i <- setNames(cbind(rownames(clustermarkers_i), clustermarkers_i, row.names = NULL), c("geneID", "myAUC", "avg_diff", "power", "avg_log2FC", "pct.1", "pct.2", "diff_exp")) 
  clustermarkers_i %>%
    top_n(n = 50, wt = avg_log2FC) ->top50
  top50<-addGeneNamesTable(top50,top50$geneID,gene_master$obgene, gene_master$merge_desc)
  
  p<-DotPlot(fil_sct_clust_findmarkers,features=top50$geneID) + RotatedAxis() + theme(axis.text.x = element_text(size = 6)) + theme(axis.text.y = element_text(size = 6)) + theme(legend.title = element_text(size = 7)) + theme(legend.text = element_text(size = 7)) + NoLegend() 
  p<-replaceLabelsX(p,gene_master$obgene,gene_master$merge_desc)
  p <- p+ coord_flip()
  p <- p+ ggtitle(paste('(avg_log2FC) top50 for cluster', m, sep = " ")) 
  print(p) 
  
}

dev.off()

