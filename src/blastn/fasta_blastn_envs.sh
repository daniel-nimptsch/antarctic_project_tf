#!/bin/bash
# .fasta_blastn_init.sh
# Custom env variables from soye
export CONDASH="/etc/profile.d/conda.sh"
export CONDABIN="/home/pipeline/.conda/envs/asv_pipeline/bin/"
export ASVBLASTDB="/home/pipeline/pipeline/localBLAST/blastdb"
export SILVADB="/home/pipeline/pipeline/localBLAST/silvadb"
export UNITEDB="/home/pipeline/pipeline/localBLAST/unite"
THREADS=24
