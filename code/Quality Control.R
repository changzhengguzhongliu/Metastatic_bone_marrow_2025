# Load R package
library(Matrix)
library(stringr)
library(tidyr)
library(magrittr)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(Seurat);packageVersion("Seurat")
library(fastSave)
# import data: liver_BM samples
dir1=list.files("/alldata",pattern = "^sample")
dir1=paste("../alldata/",dir1,sep = '')
patients=paste("LBMET",c(1:(length(dir1)/3)),sep = '')
namedir1=vector()
namedir1 <- lapply(patients, function(patient) {
  modified_patient <- paste(patient, c("Tumor", "Involved", "Distal"), sep = "-")
})
namedir1=unlist(namedir1)
names(dir1) = namedir1;dir1
# creat seurat obj
panel1 = list()
for (i in dir1) {
  panel1[[gsub("../alldata/","",i)]] =  Read10X(i)
}
seurat_obj_liver = lapply(panel1,function(x) CreateSeuratObject(counts = x,
                                                                min.cells = 3,
                                                                min.features = 100, 
                                                                project = "LIHC_bone_metastasis"))
# change orig.ident and active.ident,and add the prefix
for (i in seq_along(seurat_obj_liver)){
  seurat_obj_liver[[i]]$orig.ident = namedir1[i]
  Idents(seurat_obj_liver[[i]])= seurat_obj_liver[[i]]$orig.ident
  prefix = unique(seurat_obj_liver[[i]]$orig.ident)
  colnames(seurat_obj_liver[[i]]) = paste0(prefix,"_",colnames(seurat_obj_liver[[i]]))
}

# import data: KIRC & PRAD & BMM samples
dir2 = "/alldata"
#dir2 = list.files(dir2)[grep("^GSM[0-9]{7}_BM.*",list.files(dir2))]
dir2 = list.files(dir2)[grep("^PRAD|^BNGN|^KIRC",list.files(dir2))]
dir2=paste("../alldata/",dir2,sep = '')
panel2 <- list()
for (i in dir2) {
  b = Matrix(as.matrix(data.table::fread(i),rownames=1),sparse = T)
  panel2[[gsub("^\\.\\./alldata/|\\.count\\.csv\\.gz$","",i)]] = b
}
seurat_obj_prostate_kidney_marrow = lapply(panel2,function(x) CreateSeuratObject(counts = x,
                                                                                 min.cells = 3,
                                                                                 min.features = 100))
# Part2: Basic QC ==== 
# 2.1: calculate basic QC index
HB.genes = c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
QC_meta <- Map(function(x,y,z,w){cbind(x@meta.data,log10GenesPerUMI=y,mitoRatio=unlist(z),RBCRatio=unlist(w))}, 
               x = seurat_obj_all, 
               # number of genes detected per UMI: data set complexity
               y = lapply(seurat_obj_all, function(x){log10(x$nFeature_RNA) / log10(x$nCount_RNA)}),
               # mitochondrial ratio
               z = lapply(seurat_obj_all, function(x){PercentageFeatureSet(object = x, pattern = "^MT-") / 100}),
               # red blood cell
               w = lapply(seurat_obj_all, function(x){PercentageFeatureSet(object = x, features=HB.genes[HB.genes %in% rownames(x)]) / 100}))                                           
save(filtered_seurat, file = "../filtered_seurat.rdata")



