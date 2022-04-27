#!/bin/bash

installDIR=`grep installDIR run.configure | sed s/.*\=//`
peptide_length=`grep peptide_length run.configure | sed s/.*\=//`
dataset_name=`grep dataset_name run.configure | sed s/.*\=//`
fastq_bam_dir=`grep fastq_bam_dir run.configure | sed s/.*\=//`
output_dir=`grep output_dir run.configure | sed s/.*\=//`
list=`grep fileList run.configure | sed s/.*\=//`

cd ${installDIR}
mkdir -p "${output_dir}"
mkdir -p "${output_dir}/${dataset_name}"
logs="${output_dir}/logs"
mkdir -p "$logs"

## run HLA typing with Seq2HLA
echo "`date '+%F %T'` Pipeline started!"
./run_Seq2HLA.sh "${installDIR}" "${peptide_length}" "${dataset_name}" "${fastq_bam_dir}" "${output_dir}" >> "$logs/log.txt" 2>&1
echo "`date '+%F %T'` HLA typing has finished!"
## detect soft-clipped reads for downstream IPA analysis
./predict_polyA_spanning_reads.sh "$installDIR" "$peptide_length" "$dataset_name" "$fastq_bam_dir" "$output_dir" >> "$logs/log.txt" 2>&1
echo "`date '+%F %T'` Identification of soft-clipped reads has finished!"

## predict potential IPA neoantigens (configure the run.configure file firstly)
./run_neoepitope_pipeline.sh >> "$logs/log.txt" 2>&1

cd "${output_dir}"
for file in $(cat ${list}); do
        id=$(echo "${file}" | cut -d "." -f1)
	line_num=$(cat "${id}".reliables.fa.out|wc -l)
        if [ "${line_num}" -lt "3" ]; then
                echo ${id} has no IPA neoantigens!
        else
                echo ${id} IPA neoantigen found. Please refer to "${id}.reliables.filtered_by_uniprot.txt" file for details
        fi
done
echo "`date '+%F %T'` IPA neoantigen prediction has finished!"

cd ${output_dir}
rm headermap_* peptideSeqs* *reliables.fa* *reliables *_out* *counts.txt* *half_sequence* tmp_* key_* *_header* EncodeGencode* *ipa*table *wo *dns24 *ups24 *bed *count *summary *SAF *win48 *clusterID* *hla
