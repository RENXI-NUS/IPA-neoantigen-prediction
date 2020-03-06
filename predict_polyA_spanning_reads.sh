#!/bin/bash
current_path=$1
dataset_name=$3
bam_path=$4
path=$5
mkdir -p "${path}"/"${dataset_name}"
cd ${path}/"${dataset_name}"
list="$5/${dataset_name}.list"

for line in $(tac ${list}); do
#	if [ ! -f "${path}/${dataset_name}/${line}_potential_polyA_sites_from_softclipped_reads.txt" ]; then
		N=24
        	if (( i % N == 0 )); then
                	wait
        	fi
	        ((i++))
		(perl ${current_path}/detect_from_soft_cliping_with_bam_file_polyT_beginning_based_on_S.pl "${bam_path}/${line}" "${path}/${line}"
                perl ${current_path}/detect_from_soft_cliping_with_bam_file_polyA_end_based_on_S.pl "${bam_path}/${line}" "${path}/${line}" ) &
		#	fi
done
