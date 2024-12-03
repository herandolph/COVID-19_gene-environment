#!/bin/bash

CONDITIONS=(control control control control control control baseline_inf baseline_inf baseline_inf baseline_inf baseline_inf baseline_inf followup_inf followup_inf followup_inf followup_inf followup_inf followup_inf)
PSEUDOBULK=(pseudosums_QN)
PCS=(3 1 10 12 5 13 14 1 4 13 8 6 2 1 2 3 1 2)
TEMP_DIRS=(3expPC_control 1expPC_control 10expPC_control 12expPC_control 5expPC_control 13expPC_control 14expPC_baseline_inf 1expPC_baseline_inf 4expPC_baseline_inf 13expPC_baseline_inf 8expPC_baseline_inf 6expPC_baseline_inf 2expPC_followup_inf 1expPC_followup_inf 2expPC_followup_inf 3expPC_followup_inf 1expPC_followup_inf 2expPC_followup_inf)
OUT_DIR=(filtered_1Mb)
META=(meta_data.txt)
COUNTS=(expression_QN_by_Timepoint_noStim_noOutliers.txt)
GENE_POS=(GRCh38.92_gene_positions.txt)
CELLTYPES=(CD14_monocytes CD16_monocytes CD4_T CD8_T B NK CD14_monocytes CD16_monocytes CD4_T CD8_T B NK CD14_monocytes CD16_monocytes CD4_T CD8_T B NK)

LEN=${#CONDITIONS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
      
       CONDITION=${CONDITIONS[$NUM]}
       PC=${PCS[$NUM]}
       TEMP_DIR=${TEMP_DIRS[$NUM]}
       CELLTYPE=${CELLTYPES[$NUM]}
       
       echo " ************************************** "
       echo " BEGIN MATRIXEQTL PIPELINE FOR: $PSEUDOBULK "
       echo " CONDITION = $CONDITION;"
       echo " EXP PC TO REGRESS = $PC;"
       echo " TEMP DIRECTORY = $TEMP_DIR;"
       echo " OUTPUT DIRECTORY = $OUT_DIR;"
       echo " ************************************** "
       sbatch 4_1-2_eQTL_job_specs_1Mb.sbatch $CONDITION $PSEUDOBULK $PC $TEMP_DIR $OUT_DIR $META $COUNTS $GENE_POS $CELLTYPE &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
