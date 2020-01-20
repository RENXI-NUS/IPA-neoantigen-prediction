#!/bin/bash

file=$1
cancer=$2
cd /data/11000039/e0149673/scratch/Projects/TCGA_unsorted_bam/intron_polyadenylated_peptide_for_MS
#cd /data/11000039/e0149673/scratch/Projects/TCGA_unsorted_bam/intron_polyadenylated_peptide_for_MS/${cancer}_final_loose_criteria2_new1
for line in $(cat "${file}" | awk '{print $(NF-2)}' ); do 
#for line in $(cat "${file}" | cut -f5); do 
	echo -n ">" >> ${file}.fa
	echo "${line}" >> ${file}.fa
	echo "${line}" >> ${file}.fa
done
	
