#!/bin/bash
#
# ITSx
#
# Written by Daniel Nimptsch

fastaFile=$1
header=$2

cd ITSx/

source /home/pipeline/miniconda3/etc/profile.d/conda.sh
eval $(conda shell.bash hook)
conda activate asv_pipeline

echo /home/pipeline/miniconda3/envs/asv_pipeline/bin/ITSx -i $fastaFile -o $header --graphical T --fasta T --save-regions 'SSU,ITS1,ITS2' --cpu 8 --detailed_results T

/home/pipeline/miniconda3/envs/asv_pipeline/bin/ITSx -i $fastaFile -o $header --graphical T --fasta T --save-regions 'SSU,ITS1,ITS2' --cpu 8 --detailed_results T
