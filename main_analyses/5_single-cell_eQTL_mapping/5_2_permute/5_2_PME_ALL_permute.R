# script to run single-cell Poisson mixed effects eQTL model (permuted)

args = commandArgs(trailingOnly = TRUE)

type <- args[1]
celltype <- args[2]
cell_state <- args[3]
geneSNP_file <- args[4]
expPCs_reg <- args[5]
out_dir <- args[6]

# load libraries
library(lme4)
library(Matrix)
library(data.table)
library(plyr)
library(dplyr)

set.seed(2022)

current = getwd()
folder = "5_single-cell_eQTL_mapping"
setwd(current)

## create directory structure to save outputs
system(paste0("mkdir -p outputs/",folder,"/",type,"/",celltype,"/",out_dir,"/raw_results"))

## read in geneSNP file
geneSNPs <- read.table(paste0("inputs/",folder,"/geneSNP_files_lfsr0.1/",geneSNP_file))
genes <- geneSNPs$V1
snps <- geneSNPs$V2

## load raw UMI counts matrix and subset on correct samples
exprs_raw <- readRDS(paste0("inputs/",folder,"/expression_data/cluster_singlets/",celltype,"_raw_counts.rds"))

## load metadata 
meta <- readRDS(paste0("inputs/",folder,"/expression_data/cluster_singlets/",celltype,"_metadata.rds"))
meta$cell_ID <- rownames(meta)
meta <- subset(meta, select = c(cell_ID, nCount_RNA, batch, Dataset, percent.mt, sample_ID, individual_ID, Stim_Status, Timepoint, Age, Sex, Severity_Sampling, DSO, DUP_B2))

## filter metadata based on condition
meta$DUP_B2 <- as.character(meta$DUP_B2)
meta$DUP_B2[is.na(meta$DUP_B2)] <- "UNIQUE"
meta$DUP_B2 <- factor(meta$DUP_B2)
meta <- subset(meta, DUP_B2 %in% c("UNIQUE"))

## keep only nonstim
meta <- subset(meta, Stim_Status %in% c("noStim"))

## downsample if CD4T cells
if(celltype == "CD4_T"){
	## read in any score
	cells <- read.table(paste0("inputs/",folder,"/expression_data/exp_PCs/",celltype,"_exp_PCs_no-followups.txt"))
	exprs_raw <- exprs_raw[, colnames(exprs_raw) %in% rownames(cells)]
	meta <- meta[meta$cell_ID %in% colnames(exprs_raw),]
}
print(paste0("n cells, ", nrow(meta)[1]))

## switch genotyping IDs
meta$individual_ID <- ifelse(startsWith(meta$individual_ID, "CH"), paste0("CR",meta$individual_ID), meta$individual_ID)
meta$individual_ID <- sub('-', '', meta$individual_ID)

## subset on correct cells (only consider first timepoint COVID PATIENTS for now)
meta <- subset(meta, Timepoint %in% c("baseline_inf"))
meta$Timepoint <- factor(meta$Timepoint, levels = c("baseline_inf"))
meta$Sex <- factor(meta$Sex, levels = c("F","M"))
meta$batch <- as.factor(meta$batch)
meta$Dataset <- factor(meta$Dataset, levels = c("COVID_B1", "COVID_B2"))
meta$individual_ID <- as.factor(meta$individual_ID)
meta$sample_ID <- as.factor(meta$sample_ID)
meta$Stim_Status <- as.factor(meta$Stim_Status)

## load gene expression PCs
exp_PCs <- read.table(paste0("inputs/",folder,"/expression_data/exp_PCs/",celltype,"_exp_PCs_no-followups.txt"))
colnames(exp_PCs) <- paste0("exp",colnames(exp_PCs))
exp_PCs$cell_ID <- rownames(exp_PCs)
exp_PCs <- subset(exp_PCs, select = c("cell_ID", paste0("expPC",1:expPCs_reg)))

# load genotype dosages and geno PCs
geno <- fread(paste0("inputs/",folder,"/genotype_data/genotypes_COVID_and_IAV_merged.txt"), data.table = FALSE)
rownames(geno) <- geno$V1
geno$V1 <- NULL
## read in geno PCs
genoPCs <- read.table(paste0("inputs/",folder,"/genotype_data/genotype_PCs.txt"), header = TRUE)
colnames(genoPCs)[1] <- "individual_ID"

for(i in 1:nrow(geneSNPs)){

	gene <- genes[i]
	snp <- snps[i]

	## get snp of interest and merge with meta data
	g <- as.data.frame(t(subset(geno, rownames(geno) %in% snp)))
	g$individual_ID <- rownames(g)
	colnames(g)[1] <- "ID"
	## remove individuals with no genotyping data 
	meta_loop <- subset(meta, individual_ID %in% colnames(geno))
	meta_loop <- join(meta_loop, g, by = "individual_ID", type = "left")
	meta_loop <- subset(meta_loop, ID != "-9")

	## only keep cells in geno data in expression data (exclude follow ups)
	exprs_raw_loop <- exprs_raw[, colnames(exprs_raw) %in% meta_loop$cell_ID]

	# merge expression PCs
	meta_loop <- join(meta_loop, exp_PCs, by = "cell_ID", type = "left")

	# merge geno PCs with data
	meta_loop <- join(meta_loop, genoPCs, by = "individual_ID", type = "left")

	# load interaction covariate (e.g., batch-corrected CV from CCA)
	cov <- read.table(paste0("inputs/",folder,"/immune_scores/ssGSEA/integrate_sctransform/ssGSEA_",cell_state,"_",celltype,".txt"), header = TRUE, sep = ",")
	colnames(cov)[2] <- "cov"
	meta_loop <- join(meta_loop, cov, by = "cell_ID", type = "left")

	length(which(colnames(exprs_raw_loop)!=meta_loop$cell_ID))
	# make data frame of variables for model
	E <- as.numeric(exprs_raw_loop[gene,]) 
	G <- meta_loop$ID
	IND <- factor(meta_loop$individual_ID)
	B <- factor(meta_loop$batch)
	D <- factor(meta_loop$Dataset)
	AGE <- scale(meta_loop$Age)
	SEX <- meta_loop$Sex
	nUMI <- scale(log(meta_loop$nCount_RNA))
	MT <- meta_loop$percent.mt/100 
	## geno PCs
	subs_geno <- meta_loop[, grepl("^PC",colnames(meta_loop))]
	PC <- subs_geno[,1:3]
	## exp PCs
	expPC <- meta_loop[, grepl("expPC",colnames(meta_loop))]
	cov <- meta_loop$cov

	data <- data.frame(E,G,IND,B,D,AGE,SEX,nUMI,MT,
	                   PC1 = PC[,1], PC2 = PC[,2], PC3 = PC[,3],
	                   expPC1 = expPC[,1], expPC2 = expPC[,2], expPC3 = expPC[,3], expPC4 = expPC[,4], expPC5 = expPC[,5], cov)
	data$G <- as.numeric(as.character(data$G))
	data$cov <- as.numeric(as.character(data$cov))

	## PERMUTE IMMUNE SCORE
	data$cov_perm <- sample(data$cov)
	    
	full_model <- lme4::glmer(formula = E ~ G + (1 | IND) + (1 | B) + D + AGE + SEX + nUMI + MT + PC1 + PC2 + PC3 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + cov_perm + G*cov_perm, 
	                          family = "poisson", nAGQ = 0, data = data, control = glmerControl(optimizer = "nloptwrap"))
	null_model <- lme4::glmer(formula = E ~ G + (1 | IND) + (1 | B) + D + AGE + SEX + nUMI + MT + PC1 + PC2 + PC3 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + cov_perm, 
	                          family = "poisson", nAGQ = 0, data = data, control = glmerControl(optimizer = "nloptwrap"))
	model_lrt <- anova(null_model, full_model)

	out <- summary(full_model)$coefficients
	colnames(out) <- c("Estimate","Std.Error","zvalue","pval")
	out <- data.frame(gene=gene, snp=snp, term=row.names(out), out, lrt_pval=model_lrt$`Pr(>Chisq)`[2])
	write.table(out, paste0("outputs/",folder,"/",type,"/",celltype,"/",out_dir,"/raw_results/",gene,"_",snp,".txt"), quote = FALSE)
}

print("COMPLETED")
