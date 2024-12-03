#!/bin/bash

TYPES=(PME_apoptosis_permute PME_cholesterol_homeostasis_permute PME_coagulation_permute PME_complement_permute PME_fatty_acid_metabolism_permute PME_glycolysis_permute PME_heme_metabolism_permute PME_il2_stat5_signaling_permute PME_il6_jak_stat3_signaling_permute PME_inflammatory_response_permute PME_interferon_alpha_response_permute PME_interferon_gamma_response_permute PME_oxidative_phosphorylation_permute PME_tnfa_signaling_via_nfkb_permute)
CELLTYPE=(CD8_T)
CELLSTATES=(apoptosis cholesterol_homeostasis coagulation complement fatty_acid_metabolism glycolysis heme_metabolism il2_stat5_signaling il6_jak_stat3_signaling inflammatory_response interferon_alpha_response interferon_gamma_response oxidative_phosphorylation tnfa_signaling_via_nfkb)
GENESNP_FILE=(CD8_T_COVID_eQTL_COVID_B1_B2_IAV.txt)
PCS=(5)
OUT_DIR=(filt_5expPCs)

LEN=${#CELLSTATES[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
    	
       CELLSTATE=${CELLSTATES[$NUM]}
       TYPE=${TYPES[$NUM]}
       
       echo " BEGIN PME MODEL FOR: $CELLTYPE $CELLSTATE $GENESNP_FILE"

       sbatch 5_2_job_specs_PME_ALL_permute.sbatch $TYPE $CELLTYPE $CELLSTATE $GENESNP_FILE $PCS $OUT_DIR &
    done
wait

echo "ALL JOBS SUBMITTED AT: `date`"

## EOF ##
