library(ggplot2)
library(ggcorrplot)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(fgsea)
library(qusage)
library(data.table)
library(RColorBrewer)

current = getwd()
folder = "3_DE_analyses"
setwd(current)
system(paste0("mkdir -p outputs/",folder,"/"))
set.seed(2021)

celltypes <- c("B","CD4_T","CD8_T","CD14_monocytes","CD16_monocytes","NK")
hallmark <- read.gmt(paste0("inputs/",folder,"/genesets/h.all.v7.1.symbols.gmt"))

plot_GSEA <- function(geneset, geneset_name){

	## directories
	OUTPUTS_dir <- paste0("outputs/",folder,"/3_infection/")
	system(paste0("mkdir -p ", OUTPUTS_dir,"/fgsea_heatmaps/"))
	system(paste0("mkdir -p ", OUTPUTS_dir,"/fgsea_tables/"))
	system(paste0("mkdir -p ", OUTPUTS_dir,"/fgsea_tables/",geneset_name,"/"))

	## read in tstats
	readIn <- function(celltype){
		df <- read.table(paste0("inputs/",folder,"/infection_res/",celltype,"_results_with_qvalues.txt"), header = TRUE, sep = ",", check.names = FALSE)
		df$genes <- rownames(df)
		df <- subset(df, select = c(genes, t_stat))
		colnames(df)[2] <- paste0(celltype)
		return(df)
	}

	B <- readIn("B")
	CD4_T <- readIn("CD4_T")
	CD8_T <- readIn("CD8_T")
	NK <- readIn("NK")
	CD14_monocytes <- readIn("CD14_monocytes")
	CD16_monocytes <- readIn("CD16_monocytes")

	## combine
	t_statistics <- join(B, CD4_T, by = "genes")
	t_statistics <- join(t_statistics, CD8_T, by = "genes")
	t_statistics <- join(t_statistics, NK, by = "genes")
	t_statistics <- join(t_statistics, CD14_monocytes, by = "genes")
	t_statistics <- join(t_statistics, CD16_monocytes, by = "genes")
	t_statistics <- t_statistics[complete.cases(t_statistics), ]
	rownames(t_statistics) <-t_statistics$genes; t_statistics$genes <- NULL

	## function to calculate gsea
	for(i in 1:ncol(t_statistics)){

		name <- colnames(t_statistics)[i]

		## take the t stats of condition i and rank them
		t_i <- subset(t_statistics, select = name); colnames(t_i)[1] <- "t_stat"
		t_i$gene <- rownames(t_i)
		t_i <- data.table(t_i[,c(2,1)])
		t_i_rank <- t_i[order(t_stat), ]

		input <- setNames(t_i_rank$t_stat, t_i_rank$gene)

		fgseaRes <- fgsea(pathways = geneset, 
		                  stats = input,
		                  minSize=15,
		                  maxSize=500,
		                  nperm=100000)

		fgseaRes_csv <- data.frame(lapply(fgseaRes, as.character), stringsAsFactors = FALSE)
		write.csv(fgseaRes_csv, file = paste0(OUTPUTS_dir,"/fgsea_tables/",geneset_name,"/",name,".csv"))

		## format for downstream
		colnames(fgseaRes) <- paste0(colnames(fgseaRes),"_",name)
		fgseaRes <- fgseaRes[,c(1,3,5)]
		colnames(fgseaRes)[1] <- "pathway"
		assign(paste0(name), fgseaRes)
	}

	## combine dfs 
	merged <- merge(CD14_monocytes, CD16_monocytes, by = c("pathway"), all = FALSE)
	merged <- merge(merged, CD8_T, by = c("pathway"), all = FALSE)
	merged <- merge(merged, CD4_T, by = c("pathway"), all = FALSE)
	merged <- merge(merged, NK, by = c("pathway"), all = FALSE)
	merged <- merge(merged, B, by = c("pathway"), all = FALSE)

	## subset pathways by threshold
	thresh = 2
	merged_subset <- as.data.frame(merged[which(-log10(merged$padj_CD14_monocytes) > thresh | -log10(merged$padj_CD16_monocytes) > thresh | -log10(merged$padj_CD8_T) > thresh | -log10(merged$padj_CD4_T) > thresh | -log10(merged$padj_NK) > thresh | -log10(merged$padj_B) > thresh),])

	## melt pval
	merged_pval <- merged_subset[, grepl(paste0("padj_"), colnames(merged_subset))]; merged_pval <- cbind(merged_subset$pathway, merged_pval); colnames(merged_pval)[1] <- "pathway"
	merged_pval_m <- melt(merged_pval, id.vars = "pathway"); colnames(merged_pval_m)[2] <- "var_padj"; colnames(merged_pval_m)[3] <- "value_padj"

	## melt NES
	merged_NES <- merged_subset[, grepl(paste0("NES_"), colnames(merged_subset))]; merged_NES <- cbind(merged_subset$pathway, merged_NES); colnames(merged_NES)[1] <- "pathway"
	merged_NES_m <- melt(merged_NES, id.vars = "pathway"); colnames(merged_NES_m)[2] <- "var_NES"; colnames(merged_NES_m)[3] <- "value_NES"; merged_NES_m$pathway <- NULL

	## merge
	merged_subset_m <- cbind(merged_pval_m, merged_NES_m)

	merged_subset_m$var_NES <- gsub("CD4_T", "CD4T", merged_subset_m$var_NES)
	merged_subset_m$var_NES <- gsub("CD8_T", "CD8T", merged_subset_m$var_NES)
	merged_subset_m$pathway <- gsub("HALLMARK_", "", merged_subset_m$pathway)
	merged_subset_m <- separate(merged_subset_m, var_NES, into = c("misc","CT"), sep = "_", extra = "merge")
	merged_subset_m$CT <- factor(merged_subset_m$CT, levels = c("B","CD4T","CD8T","NK","CD14_monocytes","CD16_monocytes"))
	merged_subset_m$pathway <- gsub("_"," ",merged_subset_m$pathway)
	merged_subset_m$pathway <- tolower(merged_subset_m$pathway)

	## plot
	base_size = 4.5
	breaks_pval = c(signif(min(-log10(merged_subset_m$value_padj)), 2), signif(max(-log10(merged_subset_m$value_padj)), 2))
	limits_pval = c(signif(min(-log10(merged_subset_m$value_padj)), 2), signif(max(-log10(merged_subset_m$value_padj)), 2))
	labels_pval = c(signif(min(-log10(merged_subset_m$value_padj)), 2), signif(max(-log10(merged_subset_m$value_padj)), 2))

	breaks_NES = c(signif(min(merged_subset_m$value_NES), 1)-0.35, signif(max(merged_subset_m$value_NES), 1)+0.35)
	limits_NES = c(signif(min(merged_subset_m$value_NES), 1)-0.35, signif(max(merged_subset_m$value_NES), 1)+0.35)
	labels_NES = c(signif(min(merged_subset_m$value_NES), 1)-0.35, signif(max(merged_subset_m$value_NES), 1)+0.35)

	## color dots that are not significant gray
	merged_subset_m$value_NES_withNA <- ifelse(merged_subset_m$value_padj < 0.10, merged_subset_m$value_NES, merged_subset_m$value_NES == NA)

	plot <- ggplot(merged_subset_m, aes(x = CT, y = pathway, fill = value_NES_withNA, size = -log10(value_padj), color = CT)) + 
			geom_point(shape = 21, stroke = 0.6, color = "grey60") +
			scale_color_manual(values = rep("black", 7), name = NULL, guide = FALSE) +
	   		scale_size(name = "-log10(padj)") +
	   		#scale_fill_gradient2(guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title = "enrichment score", frame.linewidth = 1, ticks = FALSE), breaks = breaks_NES, limits = limits_NES, labels = labels_NES) +
	   		scale_fill_gradientn(colours = rev(brewer.pal(11,"Spectral")), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 0.3, title = "NES", ticks = FALSE), breaks = breaks_NES, limits = limits_NES, labels = labels_NES, na.value = "grey75") +
			guides(colour = FALSE) + 
			theme_bw(base_size = base_size)+
		  	coord_equal(ratio=0.8) +
		  	theme(legend.position = "bottom",
		  		legend.box = "vertical",
		        legend.text = element_text(size=base_size*1.5),
		        legend.title = element_text(size=base_size*1.5),
		        axis.ticks.x = element_blank(),
		        axis.ticks.y = element_blank(),
		        axis.text.x = element_text(size = base_size*1.5, angle = 45, hjust = 1, vjust = 1, colour = "black"),
		        axis.text.y = element_text(size = base_size*1.5),
		        panel.border = element_blank(), 
		        panel.grid.minor = element_blank(),
		        panel.grid.major = element_blank(),
		        #panel.grid.minor = element_line(colour="black"),
		        legend.key.height = unit(0.5, "cm"),
		        legend.key = element_rect(fill = "white"),
		        legend.key.width = unit(0.6, "cm"),
		        plot.title = element_text(size = 10)) +
		  	ggtitle(paste0(geneset_name)) + 
		  	xlab("") +
		  	ylab("")

	pdf(paste0(OUTPUTS_dir,"/fgsea_heatmaps/",geneset_name,".pdf"), width = 4, height = 7, useDingbats=FALSE)
	print(plot)
	dev.off()
}

plot_GSEA(hallmark, "hallmark")

