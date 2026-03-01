# Part1: import paper data and create seurat object ====
library(Matrix)
library(patchwork)
library(stringr) 
library(tidyr)
library(magrittr)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(Seurat,lib.loc = "/share/home/Olive/miniconda3/envs/R4.3/lib/R/library");packageVersion("Seurat")
library(harmony)
# relative parameter
resolution = seq(0.05,1.5,by=0.05);length(resolution) # The resolution to FindCluster
method_integration = c("HarmonyIntegration","RPCAIntegration",
                       "CCAIntegration","FastMNNIntegration",
                       "scVIIntegration") # The method of integration, do not change! 
new_reduction = c("harmony","integrated.rpca",
                  "integrated.cca","integrated.mnn","integrated.scvi") # The name of reduction results, not recommend to change
cluster_name_all = unlist(lapply(paste0(new_reduction,"_cluster_"),function(x){paste0(x,resolution)}))
save(cluster_name_all,file = "cluster_name_all.rdata")
DimPlot_theme = theme(axis.ticks = element_blank(),
                      axis.text = element_blank(),
                      axis.line = element_line(arrow = arrow(angle = 20, length = unit(.10,"inches"),type = "closed")))# DimPlot theme
umap_reduction_name = paste0("umap_",new_reduction) # The name of reduction results, not recommend to change
tsne_reduction_name = paste0("tsne_",new_reduction) # The name of reduction results, not recommend to change
plot_umap = NULL
plot_tsne = NULL
load("../filtered_seurat.rdata")
merge_seurat <- merge(filtered_seurat_gene_remove_doublet[[1]], filtered_seurat_gene_remove_doublet[2:length(filtered_seurat_gene_remove_doublet)])
merge_seurat <- NormalizeData(merge_seurat) %>%
                    FindVariableFeatures() %>%
                    ScaleData() %>% RunPCA()
merge_seurat <- FindNeighbors(merge_seurat, dims = 1:30, reduction = "pca")
merge_seurat <- FindClusters(merge_seurat, resolution = 2, cluster.name = "unintegrated_clusters")
merge_seurat <- RunUMAP(merge_seurat, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
merge_seurat <- RunTSNE(merge_seurat, dims = 1:30, reduction = "pca", reduction.name = "tsne.unintegrated")
for(i in 1:2){ # scVIIntegration require the conda_env of scvi-tools,so set 1:4
  print(paste0("The time is ",Sys.time(),". The integration of ",new_reduction[i]," starts")) # get the present time
  cluster_name = paste0(new_reduction[i],"_cluster_",resolution)
  merge_seurat <- IntegrateLayers(
    object = merge_seurat, method = method_integration[i],
    orig.reduction = "pca", new.reduction = new_reduction[i],verbose = FALSE) %>% 
  FindNeighbors(reduction = new_reduction[i], dims = 1:30) %>% 
  FindClusters(resolution = resolution, cluster.name = cluster_name) %>% 
  RunUMAP(reduction = new_reduction[i], dims = 1:30, reduction.name = umap_reduction_name[i]) %>% 
  RunTSNE(reduction = new_reduction[i],dims = 1:30, reduction.name = tsne_reduction_name[i])
  print(paste0("The time is ",Sys.time(),". The integration of ",new_reduction[i]," was done")) # get the present time
  save(merge_seurat,file = "merge_seurat_cycle_for.rdata")
  #load("merge_seurat_cycle_for.rdata")
  #plot
  plot_umap[[umap_reduction_name[i]]] = DimPlot(
    merge_seurat,
    reduction = umap_reduction_name[i],
    group.by = c("cancer_site", "cancer_type", "cancer_type_site",cluster_name),
    combine = FALSE, label.size = 2)
  save(plot_umap,file = "plot_umap.rdata")
  plot_tsne[[tsne_reduction_name[i]]] = DimPlot(
    merge_seurat,
    reduction = tsne_reduction_name[i],
    group.by = c("cancer_site", "cancer_type", "cancer_type_site",cluster_name),
    combine = FALSE, label.size = 2)
  save(plot_tsne,file = "plot_tsne.rdata")
  print(paste0("The time is ",Sys.time(),". The plot_umap and plot_tsne have been done(before plot unintegration)")) # get the present time
}
