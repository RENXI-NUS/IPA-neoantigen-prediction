#!/bin/bash
#PBS -N test_pipeline
#PBS -l walltime=24:00:00
#PBS -P 11000039 
#PBS -l select=1:ncpus=24

cd /data/11000039/e0149673/scratch/Projects/TCGA_unsorted_bam/scripts/test/IPA-neoantigen-prediction
curr_dir=`grep curr_dir run.config | sed s/.*\=//`
peptide_length=`grep peptide_length run.config | sed s/.*\=//`
dataset_name=`grep dataset_name run.config | sed s/.*\=//`
bam_dir=`grep bam_dir run.config | sed s/.*\=//`
output_dir=`grep output_dir run.config | sed s/.*\=//`

mkdir -p "${output_dir}"
mkdir -p "${output_dir}/${dataset_name}"
logs="${output_dir}/logs"
mkdir -p "$logs"
cd ${bam_dir}
ls *"bam" > "${output_dir}/${dataset_name}.list"

cd $curr_dir

## run HLA typing with Seq2HLA
bash "$curr_dir"/run_Seq2HLA.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out 

## detect soft-clipped reads for downstream IPA analysis
bash "$curr_dir"/predict_polyA_spanning_reads.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out 

## predict potential IPA neoantigens
bash "$curr_dir"/run_neoepitope_pipeline.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out
