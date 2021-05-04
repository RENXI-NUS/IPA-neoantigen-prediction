#!/bin/bash

cd /data/11000039/e0149673/scratch/Projects/TCGA_unsorted_bam/scripts/test/IPA-neoantigen-prediction

installDIR=`grep installDIR run.configure | sed s/.*\=//`
peptide_length=`grep peptide_length run.configure | sed s/.*\=//`
dataset_name=`grep dataset_name run.configure | sed s/.*\=//`
bam_dir=`grep bam_dir run.configure | sed s/.*\=//`
output_dir=`grep output_dir run.configure | sed s/.*\=//`

mkdir -p "${output_dir}"
mkdir -p "${output_dir}/${dataset_name}"
logs="${output_dir}/logs"
mkdir -p "$logs"
cd ${bam_dir}
ls *"bam" > "${output_dir}/${dataset_name}.list"

cd $installDIR

## run HLA typing with Seq2HLA
bash "$installDIR"/run_Seq2HLA.sh "$installDIR" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out 

## detect soft-clipped reads for downstream IPA analysis
bash "$installDIR"/predict_polyA_spanning_reads.sh "$installDIR" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out 

## predict potential IPA neoantigens (configure the run.configure file firstly)
bash "$installDIR"/run_neoepitope_pipeline.sh > log.out
