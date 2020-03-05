#!/bin/bash

set -u

bash configure.sh

pname=$(basename $0)

$result_file="$5/$3"
mkdir -p "$result_file"

[ $# -ne 5 ] && { printf "Prediction of neoantigens derived from intronic polyadenylation.
Usage: $pname [full path for the prediction] [peptide length] [cancer name] [full path for the bam file] [output path]

[current full path]		STRING without space; use _ for space.Tab file only.
[peptide length]		The length of peptide to be used for binding affinity prediction. 
[dataset name]		The dataset name for prediction. Please use TCGA abbreviation, e.g. BLCA BRCA COAD GBM LAML LIHC LUAD LUSC OV PRAD SKCM STAD
[full path for the bam files and raw fastq files]		The corresponding full path for the unsorted bam files (i.e. sort by name) and raw fastq files of the analyzed dataset. STRING without space; use _ for space.Tab file only
[outputfile path]		STRING without space; use _ for space.Tab file only.\n
"; exit 1; }

if [ $2 -lt 8 ] || [ $2 -gt 14 ]; then

printf "\nSequence length must be from 8 to 14 amino acids\n
Usage: $pname [full path for the prediction] [peptide length] [cancer name] [full path for the bam file] [output path]\n\n"

exit 1

fi

cd $1

bash "$1"/predict_polyA_spanning_reads.sh "$1" "$2" "$3" "$4" "$5" > log.out 
bash "$1"/run_Seq2HLA.sh "$1" "$2" "$3" "$4" "$5" > log.out 
bash "$1"/run_neoepitope_pipeline.sh "$1" "$2" "$3" "$4" "$5" > log.out
