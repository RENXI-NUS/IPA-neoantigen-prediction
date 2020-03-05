#!/bin/bash
current_path=$1
dataset_name=$3
bam_path=$4
path=$5
mkdir -p "${path}"/"${dataset_name}"
cd ${path}/"${dataset_name}"
list="$1/${dataset_name}.list"

for line in $(tac ${list}); do
	if [ ! -f "${path}/${dataset_name}/${line}.unsorted.bam_out_PAS_sites_from_softclipping_based_on_S.txt" ]; then
		N=24
        	if (( i % N == 0 )); then
                	wait
        	fi
	        ((i++))
		(perl ${current_path}/detect_from_soft_cliping_with_bam_file_polyT_beginning_based_on_S.pl ${bam_path}/${line}.unsorted.bam ${line}.unsorted.bam
	        perl ${current_path}/detect_from_soft_cliping_with_bam_file_polyA_end_based_on_S.pl ${bam_path}/${line}.unsorted.bam ${line}.unsorted.bam ) &
	fi
done
