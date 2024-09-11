#!/bin/bash  
#BSUB -W 72:00   # Set the wall time: HH:MM
#BSUB -n 6  ### -- ask for number of cores (default: 1) 
#BSUB -R "span[hosts=1] rusage[mem=12]"

###############################
# Arrays for different samples
###############################
IDs=(051_Primary 066_Primary 061_Primary 006_Primary)

Sample_prefixes=(ST051VFtumor ST066VFtumor ST057VFtumor ST006_FF)  

fastqs=(/beevol/home/goodspeed/core/hayashi_ewing_scRNAseq_Jul2021/select_raw/ST051VFtumor /beevol/home/goodspeed/core/hayashi_ewing_scRNAseq_Jul2021/select_raw/ST066VFtumor /beevol/home/goodspeed/core/hayashi_ewing_scRNAseq_Jul2021/select_raw/ST057VFtumor /beevol/home/goodspeed/core/hayashi_ewing_scRNAseq_Jul2021/from_bridget_sanford/10X_singlecell_data/10X-Jan2021_samples6-9-17-21/ST006_FF)
# Reading from column in Ewing_samples_8-26-21   

# Module load
module load cellranger/6.0.1

# Paths                                                                                                                                                                                                                                                               
transcriptome_10X=/beevol/home/goodspeed/genomes/human/refdata-gex-GRCh38-2020-A

# Loop through each sample
for i in "${!IDs[@]}"
do

cellranger count \
--id=${IDs[$i]} \
--jobmode=lsf \
--sample=${Sample_prefixes[$i]} \
--fastqs=${fastqs[$i]} \
--transcriptome=${transcriptome_10X}

done
