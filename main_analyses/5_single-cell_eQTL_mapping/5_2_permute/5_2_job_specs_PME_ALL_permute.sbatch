#!/bin/bash

############################################
#SBATCH --job-name=PME
#SBATCH --output=logs/PME-%j.out
#SBATCH --error=logs/PME-%j.err
#SBATCH --time=35:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-lbarreiro
#SBATCH --mem=40G
#SBATCH --nodes=1
############################################

TYPE=${1}
CELLTYPE=${2}
CELLSTATE=${3}
GENESNP_FILE=${4}
PCS=${5}
OUT_DIR=${6}

module load R/4.0.3

Rscript 5_2_PME_ALL_permute.R $TYPE $CELLTYPE $CELLSTATE $GENESNP_FILE $PCS $OUT_DIR

wait

echo " COMPLETED AT `date`"

## EOF ##
