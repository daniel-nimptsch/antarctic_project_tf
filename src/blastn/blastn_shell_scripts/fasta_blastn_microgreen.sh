#!/bin/bash
#
# Do a blastn-query with the microgreen-db and save the blastoutput to a blastoutput-file
#
# Written by Daniel Nimptsch

fastaFile=$1
blastout=$2

source ~/.fasta_blastn_envs.sh

# source ~/miniconda3/etc/profile.d/conda.sh
source $CONDASH
eval $(conda shell.bash hook)
conda activate asv_pipeline

# cd 
cd /home/soye/micro_green_db

blastn="${CONDABIN}blastn"

echo $blastn -task blastn -query $fastaFile -db microgreen_id.fasta -max_target_seqs 50 -num_threads $THREADS -outfmt "7 taxids sscinames bitscore evalue qcovs pident sacc staxids stitle" -out $blastout

$blastn -task blastn -query $fastaFile -db microgreen_id.fasta -max_target_seqs 50 -num_threads $THREADS -outfmt "7 taxids sscinames bitscore evalue qcovs pident sacc staxids stitle" -out $blastout
