library(Matrix)
library(data.table)
library(plyr)
library(dplyr)
library(qvalue)
library(mgsub)

current = getwd()
folder = "5_single-cell_eQTL_mapping"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/PME_ALL_concat_allCTs_empQ"))

in_dir <- paste0("inputs/",folder,"/PME_ALL_concat_allCTs/")
out_dir <- paste0("outputs/",folder,"/PME_ALL_concat_allCTs_empQ/")

types <- c("PME_apoptosis", "PME_cholesterol_homeostasis", "PME_coagulation", "PME_complement", "PME_fatty_acid_metabolism", "PME_glycolysis", "PME_heme_metabolism", "PME_il2_stat5_signaling", "PME_il6_jak_stat3_signaling", "PME_inflammatory_response", "PME_interferon_alpha_response", "PME_interferon_gamma_response", "PME_oxidative_phosphorylation", "PME_tnfa_signaling_via_nfkb")
types_permute <- c("PME_apoptosis_permute", "PME_cholesterol_homeostasis_permute", "PME_coagulation_permute", "PME_complement_permute", "PME_fatty_acid_metabolism_permute", "PME_glycolysis_permute", "PME_heme_metabolism_permute", "PME_il2_stat5_signaling_permute", "PME_il6_jak_stat3_signaling_permute", "PME_inflammatory_response_permute", "PME_interferon_alpha_response_permute", "PME_interferon_gamma_response_permute", "PME_oxidative_phosphorylation_permute", "PME_tnfa_signaling_via_nfkb_permute")
celltypes <- c("B","CD4_T","CD8_T","NK","CD14_monocytes","CD16_monocytes")

for(i in 1:length(types)){

	type <- types[i]
	type_permute <- types_permute[i]

	## select g * cov results only
	res <- read.table(paste0(in_dir,"concat_results_qvalue_",type,".txt"), header = TRUE)
	res_permuted <- read.table(paste0(in_dir,"concat_results_qvalue_",type_permute,".txt"), header = TRUE)

	for(j in 1:length(celltypes)){
		
		celltype_j <- celltypes[j]

		res_CT <- subset(res, celltype %in% celltype_j)
		res_perm_CT <- subset(res_permuted, celltype %in% celltype_j)

		## calculate qvalues using the permuted LRT pvalues
		emp_pvalues <- empPvals(stat = -log10(res_CT[,8]), stat0 = -log10(as.matrix(res_perm_CT[,8])), pool = T)
		qvalues <- qvalue(p = emp_pvalues)$qvalue
		res_CT <- cbind(res_CT, qvalues)
		colnames(res_CT)[11] <- "emp_qvalue"

		num <- nrow(subset(res_CT, emp_qvalue < 0.10))

		if(j == 1){
			summary_CT <- c(type, celltype_j, num)
		}else{
			summary_CT_i <- c(type, celltype_j, num)
			summary_CT <- rbind(summary_CT, summary_CT_i)
		}

		if(j == 1){
			res_final <- res_CT
		}else{
			res_final <- rbind(res_final, res_CT)
		}
	}

	write.table(res_final, file = paste0(out_dir,"/concat_results_qvalue_empP_",type,".txt"), quote = FALSE)

	if(i == 1){
		summary_PME <- summary_CT
	}else{
		summary_PME <- rbind(summary_PME, summary_CT)
	}
}

summary_PME <- as.data.frame(summary_PME)
colnames(summary_PME) <- c("model","celltype","sig_q_0.10")
write.table(summary_PME, paste0(out_dir,"PME_summary_table.txt"), quote = FALSE, row.names = FALSE)


