#!/bin/bash
#PBS -N test_pipeline
#PBS -l walltime=24:00:00
#PBS -P 11000039 
#PBS -l select=1:ncpus=24

set -u

bash configure.sh
cd $1

bash "$1"/run_Seq2HLA.sh "$1" "$2" "$3" "$4" "$5" > log.out 
bash "$1"/predict_polyA_spanning_reads.sh "$1" "$2" "$3" "$4" "$5" > log.out 
bash "$1"/run_neoepitope_pipeline.sh "$1" "$2" "$3" "$4" "$5" > log.out
