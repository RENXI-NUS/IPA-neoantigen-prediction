#!/bin/bash

set -u

pname=$(basename $0)
[ $# -ne 5 ] && { printf "Prediction of neoantigens derived from intronic polyadenylation.
Usage: $pname [full path for the prediction] [peptide length] [cancer name] [full path for the bam file] [output path]

[current full path]		STRING without space; use _ for space.Tab file only.
[peptide length]		The length of peptide to be used for binding affinity prediction.
[cancer name]		The cancer name for prediction. Please use TCGA abbreviation, e.g. BLCA BRCA COAD GBM LAML LIHC LUAD LUSC OV PRAD SKCM STAD
[full path for the bam files]		The corresponding full path for the bam files of the cancers. STRING without space; use _ for space.Tab file only
[outputfile path]		STRING without space; use _ for space.Tab file only.\n
"; exit 1; }

[ ! -f $1 ] && echo "$1 does not exist" && exit 1
[ ! -f $3 ] && echo "$3 does not exist" && exit 1
[ ! -f $4 ] && echo "$4 does not exist" && exit 1
[ -f $5 ] && echo "$5 exist. Change output file name" && exit 1

if [ $2 -lt 8 ] || [ $2 -gt 14 ]; then

printf "sequence length must be from 8 to 14 amino acids\n
Usage: $pname [full path for input fasta file] [sequence length] [number of sequence] [output file name]\n\n"

exit 1

fi

cd $1
#for cancer in "BLCA" "BRCA" "COAD" "GBM" "LAML" "LIHC" "LUAD" "LUSC" "OV" "PRAD" "SKCM" "STAD" ; do
bash "$1"/run_PASRA_TCGA.sh "$1" "$2" "$3" "$4" "$5" > log.out 
bash "$1"/run_neoepitope_pipeline.sh "$1" "$2" "$3" "$4" "$5" > log.out
#done
