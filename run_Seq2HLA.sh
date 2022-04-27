#!/bin/bash 

current_path=$1
dataset_name=$3
fastq_bam_path=$4
path=$5
python_script=`grep python2 run.configure | sed s/.*\=//`
seq2HLA=`grep seq2HLA_script run.configure | sed s/.*\=//`
list=`grep fileList run.configure | sed s/.*\=//`
cd "${path}/${dataset_name}"
rename s/fq.gz/fastq.gz/ *

for line in $(cat ${list}); do
        mkdir -p "${line}"
#	if [ ! -f "/${id}/"*"ClassI-class.HLAgenotype4digits" ]; then
#		N=2
#                if (( i % N == 0 )); then
#   	             wait
#                fi
#                ((i++))
#		(
		file1=$(find "$fastq_bam_path" -name "${line}*_1*fastq.gz")
                file2=$(find "$fastq_bam_path" -name "${line}*_2*fastq.gz")
		"$python_script" "$seq2HLA" -1 "${file1}" -2 "${file2}" -r "${path}/${dataset_name}/$line/$line" -p 80 -3 3
#		) &
#	fi
done
