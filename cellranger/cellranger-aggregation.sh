#!/bin/bash  
#BSUB -W 48:00   # Set the wall time: HH:MM
#BSUB -R "span[hosts=1] rusage[mem=12]"

# Paths
CSV=aggregation.csv

# Module load
module load cellranger/6.0.1

# Aggregation
cellranger aggr \
--id=Aggr_noNorm_12-7-21 \
--csv=${CSV} \
--normalize=none
