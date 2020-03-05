#!/bin/bash
set -u
pname=$(basename $0)
[ $# -ne 5 ] && { printf "Prediction of neoantigens derived from intronic polyadenylation.
Usage: $pname [full path for the prediction] [peptide length] [dataset name] [full path for the bam file] [output path]

[current full path]                     STRING without space; use _ for space.
[peptide length]                        The length of peptide to be used for binding affinity prediction. 
[dataset name]                          The dataset name of this prediction.
[full path to BAM and FASTQ files]      The corresponding full path for the unsorted BAM and FASTQ files. STRING without space; use _ for space.
[outputfile full path]                  STRING without space; use _ for space.\n
"; exit 1; }

if [ $2 -lt 8 ] || [ $2 -gt 14 ]; then

printf "\nSequence length must be from 8 to 14 amino acids\n
Usage: $pname [current full path] [peptide length] [dataset name] [full path to BAM and FASTQ files] [outputfile full path]\n\n"

exit 1

fi

curr_dir=$1
peptide_length=$2
dataset_name=$3
bam_dir=$4
output_dir=$5

cd $curr_dir
bash "$curr_dir"/run_Seq2HLA.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out
bash "$curr_dir"/predict_polyA_spanning_reads.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out 
bash "$curr_dir"/run_neoepitope_pipeline.sh "$curr_dir" "$peptide_length" "$dataset_name" "$bam_dir" "$output_dir" > log.out
