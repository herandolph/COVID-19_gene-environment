library(limma)
library(edgeR)
library(ggplot2)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)

current = getwd()
folder = "3_DE_analyses"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))

## inputs
celltypes <- c("B","CD4_T","CD8_T","CD14_monocytes","CD16_monocytes","NK")
residuals_dir <- paste0("inputs/",folder,"/expression_ests/")

## outputs
permutations_dir <- paste0("outputs/",folder,"/3_severity/")
system(paste0("mkdir -p ",permutations_dir))
num_permutations <- c(1:10)

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

## keep only COVID patients
md = md[which(md$Timepoint %in% c("baseline_inf")),]

## keep only nonstim
md = md[which(md$Stim_Status == "noStim"),]

## replace NA BMI with average
md$BMI <- ifelse(is.na(md$BMI), mean(md$BMI, na.rm = TRUE), md$BMI)
md$BMI_scale <- scale(md$BMI)
## replace NA comorbidities with 0
md$EC_Chronic_renal_insufficiency[is.na(md$EC_Chronic_renal_insufficiency)] <- 0
md$EC_Chronic_Heart_Insufficiency[is.na(md$EC_Chronic_Heart_Insufficiency)] <- 0
md$EC_Chronic_Lung_Insufficiency[is.na(md$EC_Chronic_Lung_Insufficiency)] <- 0
md$EC_Chronic_liver_insufficiency[is.na(md$EC_Chronic_liver_insufficiency)] <- 0
md$EC_Dyslipidemia[is.na(md$EC_Dyslipidemia)] <- 0
md$EC_DB1[is.na(md$EC_DB1)] <- 0
md$EC_DB2[is.na(md$EC_DB2)] <- 0
md$EC_Hypertension[is.na(md$EC_Hypertension)] <- 0

md$Dataset <- factor(md$Dataset, levels = c("COVID_B1", "COVID_B2"))
md$Ancestry_UMAP <- factor(md$Ancestry_UMAP)
md$EC_Chronic_renal_insufficiency <- factor(md$EC_Chronic_renal_insufficiency, levels = c("0","1"))
md$EC_Chronic_Heart_Insufficiency <- factor(md$EC_Chronic_Heart_Insufficiency, levels = c("0","1"))
md$EC_Chronic_Lung_Insufficiency <- factor(md$EC_Chronic_Lung_Insufficiency, levels = c("0","1"))
md$EC_Chronic_liver_insufficiency <- factor(md$EC_Chronic_liver_insufficiency, levels = c("0","1"))
md$EC_Hypertension <- factor(md$EC_Hypertension, levels = c("0","1"))
md$EC_Dyslipidemia <- factor(md$EC_Dyslipidemia, levels = c("0","1"))
md$EC_DB1 <- factor(md$EC_DB1, levels = c("0","1"))
md$EC_DB2 <- factor(md$EC_DB2, levels = c("0","1"))
md$Ancestry_UMAP <- factor(md$Ancestry_UMAP, levels = c("eur","afr","amr","eas","sas"))

## change severity to numeric
md$RespSupp_Sampling <- as.numeric(as.character(md$RespSupp_Sampling))

for (j in 1:length(celltypes)){

	cell_type_i <- celltypes[j]
	md_i <- md
	reorder_names <- rownames(md_i)
	
	## read in corrected expression
	residuals <- read.table(paste0(residuals_dir,cell_type_i,"_raw_expression_noStim_noOutliers.txt"), header = TRUE, sep = ",", check.names = FALSE)
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

	## add counts
	counts <- subset(md_i, select = c(paste0(cell_type_i)))
	colnames(counts)[1] <- "counts"
	md_i <- cbind(md_i, counts)

	md_i$Timepoint <- factor(md_i$Timepoint, levels = c("baseline_inf"))
	md_i$Stim_Status <- factor(md_i$Stim_Status, levels = c("noStim"))
	md_i$Dataset <- factor(md_i$Dataset, levels = c("COVID_B1","COVID_B2"))

	for (i in 1:length(num_permutations)){

		perm_number <- num_permutations[i]

		## PERMUTE SEVERITY ACROSS INDIVIDUALS
		meta_data_PERM <- md_i
		meta_data_PERM$perm_Severity <- sample(meta_data_PERM$RespSupp_Sampling)

		reorder_names <- rownames(meta_data_PERM)
		residuals <- residuals[reorder_names]
		weights <- weights[reorder_names]

		length(which(colnames(residuals)!=rownames(meta_data_PERM)))
		length(which(colnames(weights)!=rownames(meta_data_PERM)))
		length(which(colnames(residuals)!=colnames(weights)))

		## differential expression
		design = model.matrix(~ Dataset + Sex + Age_scale + BMI_scale + Ancestry_UMAP + EC_Chronic_renal_insufficiency + EC_Chronic_Heart_Insufficiency + EC_Chronic_Lung_Insufficiency + EC_Chronic_liver_insufficiency + EC_Hypertension + EC_Dyslipidemia + EC_DB1 + EC_DB2 + counts + perm_Severity, data = meta_data_PERM)
		## remove columns that are all 0s
		design <- design[, colSums(design != 0) > 0]

		vfit <- lmFit(residuals, weights = weights, design)
		vfit <- eBayes(vfit)

		pvalues = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c("perm_Severity"))]); colnames(pvalues)[1] <- paste0("pvalues_perm",perm_number,"_",cell_type_i)
		pvalues$genes <- rownames(pvalues)

		if(i == 1){
			results <- cbind(pvalues)
		}else{
			results <- merge(results, pvalues, by = "genes")
		}
	}

	write.table(results, paste0(permutations_dir,"permutations_",cell_type_i,".txt"), quote = FALSE, sep = ",", row.names = FALSE)

	print(cell_type_i)
}

