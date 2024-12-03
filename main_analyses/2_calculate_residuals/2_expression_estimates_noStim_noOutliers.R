library(limma)
library(edgeR)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)

current = getwd()
folder = "2_calculate_residuals"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## cell types
celltypes <- c("B","CD4_T","CD8_T","CD14_monocytes","CD16_monocytes","NK")

## read in meta data from individuals
md <- read.table(paste0("inputs/",folder,"/meta_data.txt"), header = TRUE)
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
rownames(md) <- md$sample_ID

## subset only on samples not stimulated 
md <- subset(md, Stim_Status == "noStim")
cell_threshold <- 5

for (i in 1:length(celltypes)){

	cell_type_i <- celltypes[i]
	meta_data_i <- md

	## remove outliers
	meta_data_i <- meta_data_i[!(meta_data_i$sample_ID %in% c("HMN171238_noStim_DSO-NA","HMN171232_noStim_DSO-NA","CHUM-1_noStim_DSO-11")),]
	meta_data_i <- meta_data_i[meta_data_i[,cell_type_i] >= cell_threshold,]

	## read in pseudobulk
	reads <- readRDS(paste0("inputs/",folder,"/",cell_type_i,"_pseudobulk.rds"))
	reorder_names <- rownames(meta_data_i)

	if(length(reorder_names) == dim(reads)[2]){
		reads <- reads[reorder_names]
	}else{
		meta_data_i <- meta_data_i[rownames(meta_data_i) %in% colnames(reads),]
		reads <- reads[,colnames(reads) %in% rownames(meta_data_i)]
		meta_data_i <- meta_data_i[rownames(meta_data_i) %in% colnames(reads),]
		reorder_names <- rownames(meta_data_i)
		reads <- reads[reorder_names]
	}

	## subset correct meta_data counts
	counts <- subset(meta_data_i, select = c(paste0(cell_type_i)))
	colnames(counts)[1] <- "counts"
	meta_data_i <- cbind(meta_data_i, counts)

	## remove lowly-expressed genes
	dge <- DGEList(counts = reads)
	dge <- calcNormFactors(dge)
	design = model.matrix(~ Dataset + Sex + Age_scale + counts, data = meta_data_i)
	## remove columns that are all 0s
	design <- design[, colSums(design != 0) > 0]
	v <- voom(dge, design, plot = TRUE)

	## median logCPM filtering
	if(cell_type_i == "CD4_T"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 1.5), ]
	}else if(cell_type_i == "CD14_monocytes" || cell_type_i == "PBMC"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 1), ]
	}else if(cell_type_i == "B" || cell_type_i == "CD8_T"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 2), ]
	}else if(cell_type_i == "NK"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 3), ]
	}else if(cell_type_i == "CD16_monocytes"){
		tab = data.frame(genes = rownames(reads), medians=apply(v$E, 1, median), order = 1:nrow(reads))
		tab = tab[order(-tab$medians), ]
		tab$order_by_median = 1:nrow(tab)
		tab = tab[order(tab$order), ]
		reads <- reads[which(tab$medians > 2.5), ]
	}

	## voom after removal of lowly expressed genes
	dge <- DGEList(counts = reads)
	dge <- calcNormFactors(dge)
	design <- model.matrix(~ Dataset + Sex + Age_scale + counts, data = meta_data_i)
	v <- voom(dge, design, plot = TRUE)

	expression <- v$E
	weights <- v$weights
	colnames(weights) <- colnames(expression)
	rownames(weights) <- rownames(expression)

	## split expression matrix into Timepoint bins (ie controls, COVID patients, recovered patients)
	quantile_norm <- function(df, cdt){

		list_CDT <- rownames(meta_data_i[meta_data_i$Timepoint %in% cdt,])
		df_QN <- df[,colnames(df) %in% list_CDT]

		## QUANTILE NORMALIZE WITHIN CONDITION
		QN_i <- matrix(, nrow = nrow(df_QN), ncol = ncol(df_QN))

		for (j in 1:nrow(QN_i)){
			exp <- df_QN[j,]
			exp_QN <- qqnorm(exp, plot = FALSE)
			exp_QN <- exp_QN$x
			QN_i[j,] <- exp_QN
		}

		rownames(QN_i) <- rownames(df_QN)
		colnames(QN_i) <- colnames(df_QN)
		return(QN_i)
	}

	QN_NI <- quantile_norm(expression, "control")
	QN_COVID <- quantile_norm(expression, "baseline_inf")
	QN_COVID_RECOV <- quantile_norm(expression, "followup_inf")
	
	length(which(rownames(QN_NI)!=rownames(QN_COVID)))
	length(which(rownames(QN_NI)!=rownames(QN_COVID_RECOV)))
	QN <- cbind(QN_NI, QN_COVID)
	QN <- cbind(QN, QN_COVID_RECOV)
	
	## write expression 
	write.table(expression, paste0("outputs/",folder,"/",cell_type_i,"_expression_noStim_noOutliers.txt"), quote = FALSE, sep = ",")
	write.table(QN, paste0("outputs/",folder,"/",cell_type_i,"_expression_QN_by_Timepoint_noStim_noOutliers.txt"), quote = FALSE, sep = ",")

	## write weights
	write.table(weights, paste0("outputs/",folder,"/",cell_type_i,"_weights_noStim_noOutliers.txt"), quote = FALSE, sep = ",")
  	print(paste0(cell_type_i, ", ", nrow(dge), " genes"))

	## correct for covariates to get GE matrix for PCA, plotting
	design <- model.matrix(~ Dataset + Sex + Age_scale + counts, data = meta_data_i)
	v <- voom(dge, design, plot = FALSE)
	vfit <-lmFit(v, design)
	vfit <- eBayes(vfit)
	corrected_expression <- v$E - vfit$coefficients[,"DatasetCOVID_B2"]%*%t(design[,"DatasetCOVID_B2"]) - vfit$coefficients[,"DatasetIAV"]%*%t(design[,"DatasetIAV"]) - vfit$coefficients[,"SexM"]%*%t(design[,"SexM"]) - vfit$coefficients[,"Age_scale"]%*%t(design[,"Age_scale"]) - vfit$coefficients[,"counts"]%*%t(design[,"counts"]) 
	write.table(corrected_expression, paste0("outputs/",folder,"/",cell_type_i,"_expression_corrected_Age_Sex_Dataset_Counts_noStim_noOutliers.txt"), quote = FALSE, sep = ",")
}

