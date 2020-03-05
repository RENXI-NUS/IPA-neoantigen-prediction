#!/bin/bash 

current_path=$1
dataset_name=$3
bam_path=$4
path=$5
mkdir -p "${path}"/"${dataset_name}"
cd ${path}/"${dataset_name}"
ls ${bam_path}/*"bam" > "${path}/${dataset_name}.list"
python_script="/data/11000039/e0149673/scratch/bin/anaconda2/envs/py2/bin/python"
seq2HLA="/data/11000039/e0149673/scratch/bin/seq2HLA/seq2HLA.py"
list="${path}/${dataset_name}.list"
logs="${path}/logs/"
mkdir -p "$logs"

for line in $(cat ${list}); do
	id=$(echo "${file}" | cut -d "." -f1)
        mkdir -p "${id}"
#	if [ ! -f "/${id}/"*"ClassI-class.HLAgenotype4digits" ]; then
		N=6
                if (( i % N == 0 )); then
                        wait
                fi
                ((i++))
		("$python_script" "$seq2HLA" -1 "${bam_path}/${id}_1.fastq.gz" -2 "${bam_path}/${id}_2.fastq.gz" -r "${path}/${dataset_name}/$id/$id" -p 12 -3 3) &
#	fi
done
