#!/bin/bash

############################################
#SBATCH --job-name=MATRIXeQTL
#SBATCH --output=logs/eQTL-%j.out
#SBATCH --error=logs/eQTL-%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=caslake
#SBATCH --mem=85G
#SBATCH --nodes=1
############################################

CONDITION=${1}
PSEUDOBULK=${2}
PC=${3}
TEMP_DIR=${4}
OUT_DIR=${5}
META=${6}
COUNTS=${7}
GENE_POS=${8}
CELLTYPE=${9}

module load R/3.6.3

Rscript 4_1-1_eQTL_mapping_100kb.R $CONDITION $PSEUDOBULK $PC $TEMP_DIR $OUT_DIR $META $COUNTS $GENE_POS $CELLTYPE

wait

echo " COMPLETED AT `date`"

## EOF ##
