# module load R/4.1.0
# R

library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(Seurat)
library(textTinyR)
library(pbapply)

current = getwd()
folder = "1_calculate_pseudobulk"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

celltypes <- c("CD4_T","CD8_T","B","CD14_monocytes","CD16_monocytes","NK")

## read in
all <- readRDS(file = paste0("inputs/",folder,"/COVID_B1_B2_IAV_Hao_et_al.rds"))
DefaultAssay(all) <- "RNA"

## split the whole object by cell type assignment because you want to do pseudobulk by cluster
all_list <- SplitObject(all, split.by = "grouped")

## loop through cell types
for (i in 1:length(all_list)){
  name <- names(all_list)[i]
  saveRDS(all_list[[i]], file = paste0("outputs/",folder,"/",name,"_cluster_singlets.rds"))
  print(name)
}

## calculate pseudobulk sums per sample for the main cell types
for (i in 1:length(celltypes)){
  name <- celltypes[i]

  ## read in single cell data (raw UMI counts stored in RNA slot) for each cluster
  dat <- readRDS(paste0("outputs/",folder,"/",name,"_cluster_singlets.rds"))
  ## get raw UMI counts
  raw_data_sparse <- GetAssayData(dat, assay = "RNA", slot = "counts")
  meta_data <- dat@meta.data
  sample_colname <- "sample_ID"

  IDs <- as.data.frame(meta_data)[, sample_colname]
  unique_ID_list <- as.list(unique(IDs))

  ## calculate pseudobulk by summing across all cells identified in each individual,cdt pair
  pseudobulk <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Sums(raw_data_sparse[,IDs == x, drop = FALSE], rowSums = TRUE)}))
  cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(raw_data_sparse[,IDs == x, drop = FALSE])})
  colnames(pseudobulk) <- names(cellcount) <- unique_ID_list
  rownames(pseudobulk) <- rownames(x = dat)

  saveRDS(pseudobulk, file = paste0("outputs/",folder,"/",name,"_pseudobulk.rds"))

  print(name)
}

