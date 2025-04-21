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
	length(which(colnames(reads)!=rownames(meta_data_i)))

	expression <- v$E
	weights <- v$weights
	colnames(weights) <- colnames(expression)
	rownames(weights) <- rownames(expression)
	write.table(expression, paste0("outputs/",folder,"/",cell_type_i,"_raw_expression_noStim_noOutliers.txt"), quote = FALSE, sep = ",")
	write.table(weights, paste0("outputs/",folder,"/",cell_type_i,"_weights_noStim_noOutliers.txt"), quote = FALSE, sep = ",")

	## split expression matrix into control vs patients so the others can be corrected for comorbidities and BMI
	CTL_samples <- subset(meta_data_i, Timepoint %in% c("control"))$sample_ID
	OTHER_samples <- subset(meta_data_i, Timepoint %in% c("baseline_inf"))$sample_ID

	exp_CTL <- expression[, colnames(expression) %in% CTL_samples]
	exp_OTHER <- expression[, colnames(expression) %in% OTHER_samples]
	weights_OTHER <- weights[, colnames(weights) %in% OTHER_samples]

	md_OTHER <- meta_data_i[rownames(meta_data_i) %in% OTHER_samples, ]
	## replace NA BMI with average
	md_OTHER$BMI <- ifelse(is.na(md_OTHER$BMI), mean(md_OTHER$BMI, na.rm = TRUE), md_OTHER$BMI)
	md_OTHER$BMI_scale <- scale(md_OTHER$BMI)
	## replace NA comorbidities with 0
	md_OTHER$EC_Chronic_renal_insufficiency[is.na(md_OTHER$EC_Chronic_renal_insufficiency)] <- 0
	md_OTHER$EC_Chronic_Heart_Insufficiency[is.na(md_OTHER$EC_Chronic_Heart_Insufficiency)] <- 0
	md_OTHER$EC_Chronic_Lung_Insufficiency[is.na(md_OTHER$EC_Chronic_Lung_Insufficiency)] <- 0
	md_OTHER$EC_Chronic_liver_insufficiency[is.na(md_OTHER$EC_Chronic_liver_insufficiency)] <- 0
	md_OTHER$EC_Dyslipidemia[is.na(md_OTHER$EC_Dyslipidemia)] <- 0
	md_OTHER$EC_DB1[is.na(md_OTHER$EC_DB1)] <- 0
	md_OTHER$EC_DB2[is.na(md_OTHER$EC_DB2)] <- 0
	md_OTHER$EC_Hypertension[is.na(md_OTHER$EC_Hypertension)] <- 0

	md_OTHER$Dataset <- factor(md_OTHER$Dataset, levels = c("COVID_B1", "COVID_B2"))
	md_OTHER$Ancestry_UMAP <- factor(md_OTHER$Ancestry_UMAP)
	md_OTHER$EC_Chronic_renal_insufficiency <- factor(md_OTHER$EC_Chronic_renal_insufficiency, levels = c("0","1"))
	md_OTHER$EC_Chronic_Heart_Insufficiency <- factor(md_OTHER$EC_Chronic_Heart_Insufficiency, levels = c("0","1"))
	md_OTHER$EC_Chronic_Lung_Insufficiency <- factor(md_OTHER$EC_Chronic_Lung_Insufficiency, levels = c("0","1"))
	md_OTHER$EC_Chronic_liver_insufficiency <- factor(md_OTHER$EC_Chronic_liver_insufficiency, levels = c("0","1"))
	md_OTHER$EC_Hypertension <- factor(md_OTHER$EC_Hypertension, levels = c("0","1"))
	md_OTHER$EC_Dyslipidemia <- factor(md_OTHER$EC_Dyslipidemia, levels = c("0","1"))
	md_OTHER$EC_DB1 <- factor(md_OTHER$EC_DB1, levels = c("0","1"))
	md_OTHER$EC_DB2 <- factor(md_OTHER$EC_DB2, levels = c("0","1"))
	md_OTHER$Ancestry_UMAP <- factor(md_OTHER$Ancestry_UMAP, levels = c("eur","afr","amr","eas","sas"))

	print(length(which(colnames(exp_OTHER)!=rownames(md_OTHER))))
	print(length(which(colnames(weights_OTHER)!=rownames(md_OTHER))))


	## correct exp_OTHER for comorbidities, BMI 
	design2 <- model.matrix(~ BMI_scale + EC_Chronic_renal_insufficiency + EC_Chronic_Heart_Insufficiency + EC_Chronic_Lung_Insufficiency + EC_Chronic_liver_insufficiency + EC_Hypertension + EC_Dyslipidemia + EC_DB1 + EC_DB2, data = md_OTHER)
	vfit <- lmFit(exp_OTHER, weights = weights_OTHER, design2)
	vfit <- eBayes(vfit)
	residuals2 <- residuals.MArrayLM(object = vfit, exp_OTHER)
	intercept <- vfit$coefficients[,1]
	corrected_expression_OTHER <- apply(residuals2, 2, function(x){x + intercept}) 

	## combine with non-corrected exp_CTL to model infection effects later
	print(length(which(rownames(corrected_expression_OTHER)!=rownames(exp_CTL))))
	combined <- cbind(corrected_expression_OTHER, exp_CTL)
	write.table(combined, paste0("outputs/",folder,"/",cell_type_i,"_expression_uncorrected-CTL_corrected_BMI-CoMorbid-PATIENTS.txt"), quote = FALSE, sep = ",")
}

