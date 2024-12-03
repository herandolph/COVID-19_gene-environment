#!/bin/bash

TYPES=(PME_apoptosis PME_cholesterol_homeostasis PME_coagulation PME_complement PME_fatty_acid_metabolism PME_glycolysis PME_heme_metabolism PME_il2_stat5_signaling PME_il6_jak_stat3_signaling PME_inflammatory_response PME_interferon_alpha_response PME_interferon_gamma_response PME_oxidative_phosphorylation PME_tnfa_signaling_via_nfkb)
CELLTYPE=(CD14_monocytes)
CELLSTATES=(apoptosis cholesterol_homeostasis coagulation complement fatty_acid_metabolism glycolysis heme_metabolism il2_stat5_signaling il6_jak_stat3_signaling inflammatory_response interferon_alpha_response interferon_gamma_response oxidative_phosphorylation tnfa_signaling_via_nfkb)
GENESNP_FILE=(CD14_monocytes_COVID_eQTL_COVID_B1_B2_IAV.txt)
PCS=(5)
OUT_DIR=(filt_5expPCs)

LEN=${#CELLSTATES[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
    	
       CELLSTATE=${CELLSTATES[$NUM]}
       TYPE=${TYPES[$NUM]}
       
       echo " BEGIN PME MODEL FOR: $CELLTYPE $CELLSTATE $GENESNP_FILE"

       sbatch 5_1_job_specs_PME_ALL.sbatch $TYPE $CELLTYPE $CELLSTATE $GENESNP_FILE $PCS $OUT_DIR &
    done
wait

echo "ALL JOBS SUBMITTED AT: `date`"

## EOF ##
