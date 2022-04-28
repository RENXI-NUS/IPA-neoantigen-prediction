#!/bin/bash
current_path=$1
dataset_name=$3
fastq_bam_path=$4
path=$5
STAR=`grep Star_tool run.configure | sed s/.*\=//`
STAR_library=`grep Star_Lib run.configure | sed s/.*\=//`
list=`grep fileList run.configure | sed s/.*\=//`
mkdir -p "${path}"/"${dataset_name}"
cd ${path}/"${dataset_name}"

for line in $(tac ${list}); do
#	if [ ! -f "${path}/${dataset_name}/${line}_potential_polyA_sites_from_softclipped_reads.txt" ]; then
#		N=2
#        	if (( i % N == 0 )); then
#                	wait
#        	fi
#	        ((i++))
#		(
                file1=$(find "$fastq_bam_path" -type f -name "${line}*_1*fastq.gz")
                file2=$(find "$fastq_bam_path" -type f -name "${line}*_2*fastq.gz")
		"$STAR" --genomeDir "$STAR_library" --runThreadN 20 --limitBAMsortRAM 100000000000 --limitIObufferSize 1000000000 1000000000 --readFilesIn "$file1" "$file2" --outSAMtype BAM Unsorted --outFileNamePrefix "$path/$line" --readFilesCommand zcat
		perl ${current_path}/detect_from_soft_cliping_with_bam_file_polyT_beginning_based_on_S.pl "${path}/${line}Aligned.out.bam" "${path}/${dataset_name}/${line}"
                perl ${current_path}/detect_from_soft_cliping_with_bam_file_polyA_end_based_on_S.pl "${path}/${line}Aligned.out.bam" "${path}/${dataset_name}/${line}" 
#		) &
#		fi
done
