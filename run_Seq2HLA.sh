#!/bin/bash 

current_path=$1
dataset_name=$3
bam_path=$4
path=$5
cd ${current_path} 
python_script=`grep py2 run.configure | sed s/.*\=//`
seq2HLA=`grep seq2HLA_script run.configure | sed s/.*\=//`

cd "${path}/${dataset_name}"
list="${path}/${dataset_name}.list"

for line in $(cat ${list}); do
	id=$(echo "${line}" | cut -d "." -f1)
        mkdir -p "${id}"
#	if [ ! -f "/${id}/"*"ClassI-class.HLAgenotype4digits" ]; then
#		N=6
#               if (( i % N == 0 )); then
#                        wait
#               fi
#               ((i++))
#		(
		"$python_script" "$seq2HLA" -1 "${bam_path}/${id}_1.fastq.gz" -2 "${bam_path}/${id}_2.fastq.gz" -r "${path}/${dataset_name}/$id/$id" -p 12 -3 3
#		) &
#	fi
done
