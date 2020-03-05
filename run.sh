#!/bin/bash
set -u

bash configure.sh
cd $curr_dir

bash "$1"/run_Seq2HLA.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out
bash "$1"/predict_polyA_spanning_reads.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out
bash "$1"/run_neoepitope_pipeline.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out
