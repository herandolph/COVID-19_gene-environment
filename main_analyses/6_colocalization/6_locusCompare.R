library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(locuscomparer)
library(stringi)
library(stringr)
library(data.table)

current = getwd()
folder = "6_colocalization"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/coloc_plots"))

windowSize <- 1000e3
geneSNPs <- c("IFNAR2_21_33237639_A_G","SNRPD2_19_45687273_A_G")
celltypes <- c("CD4_T","CD14_monocytes")
GWAS <- c("B2","B2")
conditions <- c("baseline_inf","baseline_inf")
out_dir <- paste0("outputs/",folder,"/coloc_plots/")

for(i in 1:length(geneSNPs)){

	celltype_i <- celltypes[i]
	geneSNP_i <- geneSNPs[i]
	GWAS_i <- GWAS[i]

	gwas <- as.data.frame(fread(paste0("inputs/",folder,"/GWAS/",GWAS_i,"_EUR.tsv.gz")))
	colnames(gwas)[1] <- c("chr")
	rsids <- as.data.frame(fread(paste0("inputs/",folder,"/rsids_chrPos_GRCh38_noHeader.txt")))
	colnames(rsids) <- c("chr","position","rsid")
	rsids$chr_pos <- paste0(rsids$chr,"_",rsids$position)
	rsids$chr <- NULL; rsids$position <- NULL

	for(j in 1:length(conditions)){

		condition_j <- conditions[j]

		## read in eqtl results
		eqtl <- as.data.frame(fread(Sys.glob(paste0("inputs/",folder,"/eQTL_results/result_original_PCs*_",celltype_i,"_",condition_j,"_1Mb.txt"))))

		tmpLoci <- geneSNP_i ## top SNP
		tmpGene <- strsplit(geneSNP_i, "_")[[1]][1]
		tmpChr <- as.integer(strsplit(geneSNP_i, "_")[[1]][2])
		tmpPos <- as.integer(strsplit(geneSNP_i, "_")[[1]][3])
		upStream <- tmpPos + windowSize
		dnStream <- tmpPos - windowSize

		## subset eqtl
		tmpEqtl <- eqtl[grep(paste0("^",tmpChr,"_"), eqtl$snps),]
		tmpEqtl <- subset(tmpEqtl, gene %in% tmpGene)
		tmpEqtl <- separate(tmpEqtl, snps, "_", into = c("chr","pos_BP"), extra = "merge", remove = FALSE)
		tmpEqtl <- separate(tmpEqtl, pos_BP, "_", into = c("position","BP"), extra = "merge", remove = TRUE)
		tmpEqtl <- tmpEqtl %>% filter(chr == tmpChr & position < upStream & position > dnStream)
		tmpEqtl$chr_pos <- paste0(tmpEqtl$chr,"_",tmpEqtl$position)

		## add rsid
		tmpEqtl_rsid <- join(tmpEqtl, rsids, by = "chr_pos", type = "left")
		## define top snp
		top_rsid <- subset(tmpEqtl_rsid, chr_pos %in% paste0(tmpChr,"_",tmpPos))$rsid[1]
		tmpEqtl_rsid <- subset(tmpEqtl_rsid, select = c(rsid, pval = pvalue))
		colnames(tmpEqtl_rsid)[2] <- "pval"
		print(paste0("n eqtl genes: ",nrow(tmpEqtl_rsid)))

		## subset gwas
		tmpGwas <- gwas %>% filter(chr == tmpChr & POS < upStream & POS > dnStream)
		tmpGwas <- subset(tmpGwas, select = c(rsid, pval = all_inv_var_meta_p))
		colnames(tmpGwas)[2] <- "pval"

		system(paste0("mkdir -p ",out_dir,"inputs/",tmpGene,"/"))
		write.table(tmpEqtl_rsid, paste0(out_dir,"inputs/",tmpGene,"/eqtl_",condition_j,".tsv"), sep = "\t", quote = FALSE)
		write.table(tmpGwas, paste0(out_dir,"inputs/",tmpGene,"/gwas_",condition_j,".tsv"), sep = "\t", quote = FALSE)

		eqtl_path <- paste0(out_dir,"inputs/",tmpGene,"/eqtl_",condition_j,".tsv")
		gwas_path <- paste0(out_dir,"inputs/",tmpGene,"/gwas_",condition_j,".tsv")

		pdf(paste0(out_dir,celltype_i,"_",tmpGene,"_",condition_j,".pdf"), height = 3.5, width = 6.5)
		print(locuscompare(in_fn1 = gwas_path, in_fn2 = eqtl_path, title = 'GWAS', title2 = 'eQTL', snp = top_rsid, genome = "hg38", population = "EUR"))
		dev.off()

		print(paste0("locus: ", tmpLoci,", condition: ",condition_j,", n eqtl genes: ",nrow(tmpEqtl_rsid),", n GWAS: ", nrow(gwas_path)))
	}
}

