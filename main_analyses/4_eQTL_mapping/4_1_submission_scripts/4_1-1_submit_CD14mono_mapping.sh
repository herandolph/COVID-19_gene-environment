#!/bin/bash

CONDITIONS=(control control control control control control control control control control baseline_inf baseline_inf baseline_inf baseline_inf baseline_inf baseline_inf baseline_inf baseline_inf baseline_inf baseline_inf)
PSEUDOBULK=(pseudosums_QN)
PCS=(1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10)
TEMP_DIRS=(1expPC_control 2expPC_control 3expPC_control 4expPC_control 5expPC_control 6expPC_control 7expPC_control 8expPC_control 9expPC_control 10expPC_control 1expPC_baseline_inf 2expPC_baseline_inf 3expPC_baseline_inf 4expPC_baseline_inf 5expPC_baseline_inf 6expPC_baseline_inf 7expPC_baseline_inf 8expPC_baseline_inf 9expPC_baseline_inf 10expPC_baseline_inf)
OUT_DIRS=(filtered_1expPCs filtered_2expPCs filtered_3expPCs filtered_4expPCs filtered_5expPCs filtered_6expPCs filtered_7expPCs filtered_8expPCs filtered_9expPCs filtered_10expPCs filtered_1expPCs filtered_2expPCs filtered_3expPCs filtered_4expPCs filtered_5expPCs filtered_6expPCs filtered_7expPCs filtered_8expPCs filtered_9expPCs filtered_10expPCs)
META=(meta_data.txt)
COUNTS=(expression_QN_by_Timepoint_noStim_noOutliers.txt)
GENE_POS=(GRCh38.92_gene_positions.txt)
CELLTYPE=(CD14_monocytes)

LEN=${#CONDITIONS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
    	
       CONDITION=${CONDITIONS[$NUM]}
       PC=${PCS[$NUM]}
       TEMP_DIR=${TEMP_DIRS[$NUM]}
       OUT_DIR=${OUT_DIRS[$NUM]}
       

       echo " ************************************** "
       echo " BEGIN MATRIXEQTL PIPELINE FOR: $PSEUDOBULK "
       echo " CONDITION = $CONDITION;"
       echo " EXP PC TO REGRESS = $PC;"
       echo " TEMP DIRECTORY = $TEMP_DIR;"
       echo " OUTPUT DIRECTORY = $OUT_DIR;"
       echo " ************************************** "
       sbatch 4_1-1_eQTL_job_specs_100kb.sbatch $CONDITION $PSEUDOBULK $PC $TEMP_DIR $OUT_DIR $META $COUNTS $GENE_POS $CELLTYPE &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
