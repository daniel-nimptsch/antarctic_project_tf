#!/bin/bash
#
# Do a blastn-query with the nt-db and save the blastoutput to a blastoutput-file
#
# Written by Daniel Nimptsch

fastaFile=$1
blastout=$2

source ~/.fasta_blastn_init.sh

echo $CONDASH
echo $ASVBLASTDB
echo $THREADS 

# source ~/miniconda3/etc/profile.d/conda.sh
source $CONDASH
eval $(conda shell.bash hook)
conda activate asv_pipeline

# cd /home/pipeline/bigdata/localBLAST/blastdb/
cd $ASVBLASTDB

echo ~/miniconda3/envs/asv_pipeline/bin/blastn -task blastn -query $fastaFile -db nt -max_target_seqs 100 -num_threads $THREADS -outfmt "7 taxids sscinames bitscore evalue qcovs pident sacc staxids stitle" -out $blastout

~/miniconda3/envs/asv_pipeline/bin/blastn -task blastn -query $fastaFile -db nt -max_target_seqs 100 -num_threads $THREADS -outfmt "7 taxids sscinames bitscore evalue qcovs pident sacc staxids stitle" -out $blastout
