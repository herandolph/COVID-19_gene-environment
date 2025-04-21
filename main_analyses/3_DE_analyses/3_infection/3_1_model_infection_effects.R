# module load R/4.1.0
# R

library(limma)
library(edgeR)
library(ggplot2)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(qvalue)
library(cobs)

current = getwd()
folder = "3_DE_analyses"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

setwd("common_functions")
source("permFDR.R")
setwd(current)

## inputs
celltypes <- c("B","CD4_T","CD8_T","CD14_monocytes","CD16_monocytes","NK")
residuals_dir <- paste0("inputs/",folder,"/expression_ests/")
permutations_dir <- paste0("inputs/",folder,"/infection_permutations/")

md <- read.table(paste0("inputs/",folder,"/meta_data.txt"), header = TRUE, sep = ",")
md$sample_ID <- as.character(md$sample_ID)
md$individual_ID <- factor(md$individual_ID)
md$Sex <- factor(md$Sex, levels = c("F","M"))
md$Severity_Sampling <- factor(md$Severity_Sampling, levels = c("CTL","MODERATE","SEVERE","CRITICAL","FOLLOWUP"))
md$RespSupp_Sampling <- factor(md$RespSupp_Sampling, levels = c("0","1","2","3","4","-1"))
md$Dataset <- factor(md$Dataset, levels = c("COVID_B1", "COVID_B2", "IAV"))
md$Timepoint <- factor(md$Timepoint, levels = c("control", "baseline_inf", "followup_inf"))
md$Stim_Status <- factor(md$Stim_Status, levels = c("noStim", "LPS"))
md$Age_scale <- scale(md$Age)
md$DUP_B2 <- factor(md$DUP_B2)
md$Ancestry_UMAP <- factor(md$Ancestry_UMAP, levels = c("eur","afr","amr","eas","sas"))
rownames(md) <- md$sample_ID

## remove duplicated
md$DUP_B2 <- as.character(md$DUP_B2)
md$DUP_B2[is.na(md$DUP_B2)] <- "UNIQUE"
md$DUP_B2 <- factor(md$DUP_B2)
md = md[which(md$DUP_B2 == "UNIQUE"),]

## keep only certain condition
md = md[which(md$Timepoint %in% c("control", "baseline_inf")),]

## keep only nonstim
md = md[which(md$Stim_Status == "noStim"),]

for (i in 1:length(celltypes)){

	cell_type_i <- celltypes[i]
	md_i <- md
	reorder_names <- rownames(md_i)

	## read in corrected expression
	residuals <- read.table(paste0(residuals_dir,cell_type_i,"_expression_uncorrected-CTL_corrected_BMI-CoMorbid-PATIENTS.txt"), header = TRUE, sep = ",", check.names = FALSE)
	residuals <- residuals[colnames(residuals) %in% md_i$sample_ID]

	## read in weights
	weights <- read.table(paste0(residuals_dir,cell_type_i,"_weights_noStim_noOutliers.txt"), header = TRUE, sep = ",", check.names = FALSE)
	weights <- weights[colnames(weights) %in% md_i$sample_ID]

	if(length(reorder_names) == dim(residuals)[2]){
		residuals <- residuals[reorder_names]
		weights <- weights[reorder_names]
	}else{
		correct_names <- colnames(residuals)
		md_i <- md_i[rownames(md_i) %in% correct_names,]
		weights <- weights[,colnames(weights) %in% correct_names]
		reorder_names <- rownames(md_i)
		residuals <- residuals[reorder_names]
		weights <- weights[reorder_names]
	}

	# read in permutations
	permutations <- read.table(paste0(permutations_dir,"permutations_",cell_type_i,".txt"), header = TRUE, sep = ",")

	## add counts
	counts <- subset(md_i, select = c(paste0(cell_type_i)))
	colnames(counts)[1] <- "counts"
	md_i <- cbind(md_i, counts)

	md_i$Timepoint <- factor(md_i$Timepoint, levels = c("control", "baseline_inf"))
	md_i$Stim_Status <- factor(md_i$Stim_Status, levels = c("noStim"))

	length(which(rownames(md_i)!=colnames(residuals)))
	length(which(colnames(weights)!=colnames(residuals)))

	## model infection differential expression
	design = model.matrix(~ Dataset + Sex + Age_scale + Ancestry_UMAP + counts + Timepoint, data = md_i)
	design <- design[, colSums(design != 0) > 0]

	vfit <- lmFit(residuals, weights = weights, design)
	vfit <- eBayes(vfit)

	## collect outputs
	betas = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("Timepointbaseline_inf"))]); colnames(betas)[1] <- "betas"
	p_values = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c("Timepointbaseline_inf"))]); colnames(p_values)[1] <- "pvalues"
	fdrs = as.data.frame(p.adjust(p_values[,1], method = "BH")); colnames(fdrs)[1] <- "fdrs"
	t_stats = as.data.frame(cbind(rownames(topTable(vfit, coef = "Timepointbaseline_inf", number = Inf)), topTable(vfit, coef = "Timepointbaseline_inf", number = Inf)$t))
	colnames(t_stats) <- c("genes","t_stat")

	results <- cbind(betas, p_values, fdrs)
	results$genes <- rownames(results)
	results <- join(results, t_stats, by = "genes")

	## match row names for permutations
	perms <- permutations[match(results$genes, permutations$genes),]
	length(which(results$genes!=perms$genes))
	rownames(results) <- results$genes; results$genes <- NULL
	rownames(perms) <- perms$genes; perms$genes <- NULL
	
	## calculate qvalues using a permutation-based null
	perm_fdrs = permFDR(full_data = results, full_column_id = "pvalues", perm_data = perms, perm_column_ids = "all", output_name = paste0("outputs/",folder,"/3_infection/"))
	write.table(perm_fdrs$fdrs, file = paste0("outputs/",folder,"/3_infection/",cell_type_i,"_results_with_qvalues.txt"), quote = FALSE, sep = ",")

	corr_fdrs <- subset(perm_fdrs$fdrs)

  	if(i == 1){
  		numGenes <- dim(results)[1]
  		numDEgenes_10 <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.10 & abs(corr_fdrs[,"betas"]) > 0.5))
  		numDEgenes_05 <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.05 & abs(corr_fdrs[,"betas"]) > 0.5))
  		numDEgenes_05_no_cutoff <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.05))
  		numDEgenes_01 <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.01 & abs(corr_fdrs[,"betas"]) > 0.5))
  		celltype <- cell_type_i
  	}else{
  		numGenes[i] <- dim(results)[1]
  		numDEgenes_10[i] <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.10 & abs(corr_fdrs[,"betas"]) > 0.5))
  		numDEgenes_05[i] <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.05 & abs(corr_fdrs[,"betas"]) > 0.5))
  		numDEgenes_05_no_cutoff[i] <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.05))
  		numDEgenes_01[i] <- length(which(corr_fdrs[,"Fdr_SAB"] < 0.01 & abs(corr_fdrs[,"betas"]) > 0.5))
  		celltype[i] <- cell_type_i
  	}

  	print(cell_type_i)
}

numDEgenes_per_CT <- as.data.frame(cbind(celltype, numDEgenes_10, numDEgenes_05, numDEgenes_01, numDEgenes_05_no_cutoff, numGenes))
colnames(numDEgenes_per_CT) <- c("celltype","numDEgenes_fdr10","numDEgenes_fdr05","numDEgenes_fdr01","numDEgenes_fdr05_no_logFC_cutoff","num_genes_tested")
write.table(numDEgenes_per_CT, paste0("outputs/",folder,"/3_infection/number_of_DEG_by_CT_qvalues_logFCcutoff.txt"), quote = FALSE)
