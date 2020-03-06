#!/bin/bash

path=$1
window=$2
cancer=$3
bam_path=$4
file_path=$5

cd ${path}
cp findPeak.cluster.pl splicing_coordinates_TCGA_12_cancers1 EncodeGencode_merged_intron_selected_excluded_exon_- EncodeGencode_merged_intron_selected_excluded_exon_+ callCluster.pl gencode_PASs_extended_by_100_filtered hg19_intron.excluded.exon.bed PASRA2PeptideSeqs.py ${file_path}
cd ${file_path}
python=/data/11000039/e0149673/scratch/bin/anaconda2/envs/py2/bin/python
netMHCpan=/data/11000039/e0149673/scratch/bin/software/netMHCpan-4.0/netMHCpan

## identify intronic PASs 
files="${cancer}.list"
for file in $(cat ${files}); do awk '$4 !~ /N/' ${file}_potential_polyA_sites_from_softclipped_reads.txt >> ${file}.0; done ##remove splicing reads
for file in $(cat ${files}); do awk 'NR==FNR{a[$2]++; next} a[$2]>=2' <(sort -k3,3r -k1,1 -k2,2n ${file}.0) <(sort -k3,3r -k1,1 -k2,2n ${file}.0) > ${file}.0.1; done	# print the PASs that supported by at least two soft-clipping reads
for file in $(cat ${files}); do sort -k3,3r -k1,1 -k2,2n ${file}.0.1 | awk '{if ($3 == "+") print $1"\t"$2"\t"$2+1"\t"$4"\t0\t"$3}' > ${file}.0.1.output; done
for file in $(cat ${files}); do sort -k3,3r -k1,1 -k2,2n ${file}.0.1 | awk '{if ($3 == "-") print $1"\t"$2-1"\t"$2"\t"$4"\t0\t"$3}' >> ${file}.0.1.output; done
for file in $(cat ${files}); do perl callCluster.pl ${file}.0.1.output ${file}.0.1.output.cluster ${file}.0.1.output.length ; done
for file in $(cat ${files}); do awk '{print $1"\t"$2"\t"$3"\t"$4":"FILENAME"\t"$5"\t"$6}' ${file}.0.1.output.cluster | sed 's/\.unsorted.bam_out_PAS_sites_from_softclipping_based_on_S.txt.0.1.output.cluster//g' | awk '{ if ($1 == "1" || $1 == "2" || $1 == "3" || $1 == "4" ||$1 =="5"||$1 == "6" ||$1 =="7" || $1 =="8"||$1=="9"||$1=="10"||$1=="11"||$1=="12"||$1=="13"||$1=="14"||$1=="15"||$1=="16"||$1=="17"||$1=="18"||$1=="19"||$1=="20"||$1=="21"||$1=="22"||$1=="X"||$1=="Y") print}' >> ${cancer}.cluster.bed ; done 
sort -k1,1 -k2,2n ${cancer}.cluster.bed > ${cancer}.cluster.sorted.bed
mergeBed -d 24 -s -i ${cancer}.cluster.sorted.bed -c 4,5,6 -o collapse,sum,distinct > ${cancer}.cluster.merged.bed
for file in $(cat ${files}); do cat ${file}.0.1.output >> ${cancer}_all_SC_reads ; done
perl findPeak.cluster.pl ${cancer}.cluster.merged.bed ${cancer}_all_SC_reads ${cancer}.cluster.merged_withPeak.bed
awk -F ',' '{if (NF >= 1 ) print }' ${cancer}.cluster.merged_withPeak.bed | awk '{print "chr"$1"\t"$7"\t"$8"\t"$4"\t"$5"\t"$6}' > ${cancer}.cluster.merged_withPeak.reliable.bed ##only keep the polyA clusters that supported by at least 2 samples. But if there is only one sample to be analyzed, then change it to one
bedtools intersect -wa -wb -a ${cancer}.cluster.merged_withPeak.reliable.bed -b B.hg19_intron.excluded.exon.bed | awk '{if (($6 == "+" && $12 == "+") || ($6 == "-" && $12 == "-")) print}' |awk -F"\t" '!seen[$1, $2, $3, $6]++' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ${cancer}.cluster.merged_withPeak.reliable.bed_only_in_intron
bedtools intersect -wa -wb -v -a <(sort -k6,6r -k1,1 -k2,2n ${cancer}.cluster.merged_withPeak.reliable.bed_only_in_intron) -b <(sort -k6,6r -k1,1 -k2,2n gencode_PASs_extended_by_100_filtered ) > ${cancer}.cluster.merged_withPeak.reliable.bed_only_in_intron_filtered1  ## filter the identified IPA with GENCODE annotation
bedtools intersect -wa -wb -v -a <(sort -k6,6r -k1,1 -k2,2n ${cancer}.cluster.merged_withPeak.reliable.bed_only_in_intron_filtered1) -b <(sort -k6,6r -k1,1 -k2,2n A.splicing_coordinates_TCGA_12_cancers ) > ${cancer}.cluster.merged_withPeak.reliable.bed_only_in_intron_filtered3  ## filter out the coordinates if has at least 6 junction reads
bedtools intersect -wa -wb -a <(sort -k6,6r -k1,1 -k2,2n ${cancer}.cluster.merged_withPeak.reliable.bed_only_in_intron_filtered3) -b <(sort -k6,6r -k1,1 -k2,2n C.EncodeGencode_merged_intron_selected_excluded_exon_+) | awk '{if (($6 == "+" && $12 == "+") || ($6 == "-" && $12 == "-")) print}' |awk -F"\t" '!seen[$1, $2, $3, $6]++' > EncodeGencode_merged_intron_selected_+_overlapped_with_${cancer}
bedtools intersect -wa -wb -a <(sort -k6,6r -k1,1 -k2,2n ${cancer}.cluster.merged_withPeak.reliable.bed_only_in_intron_filtered3) -b <(sort -k6,6r -k1,1 -k2,2n D.EncodeGencode_merged_intron_selected_excluded_exon_-) | awk '{if (($6 == "+" && $12 == "+") || ($6 == "-" && $12 == "-")) print}' |awk -F"\t" '!seen[$1, $2, $3, $6]++' > EncodeGencode_merged_intron_selected_-_overlapped_with_${cancer}
awk -F'\t' -v OFS='\t' '{ if ($6 == "+") print $1,$2,$10,$13,$14,$15,$6,$8,$5,$11}' EncodeGencode_merged_intron_selected_+_overlapped_with_${cancer} > input_for_peptideseqs_of_1_${cancer}
awk -F'\t' -v OFS='\t' '{ if ($6 == "-") print $1,$3,$10,$13,$14,$15,$6,$9,$5,$11}' EncodeGencode_merged_intron_selected_-_overlapped_with_${cancer} >> input_for_peptideseqs_of_1_${cancer}
awk -F"\t" '!seen[$1, $2, $(NF-2)]++' input_for_peptideseqs_of_1_${cancer} > ${cancer}_input_for_peptideseqs_uniq

## generate IPA derived peptide sequences
${python} GeneratePeptide.py "${window}" "${file_path}" ${cancer}_input_for_peptideseqs_uniq
grep '>s=' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq.fa > peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header

## filter out the IPA without significant coverage drop between upstream and downstream of the IPA site
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,"=|:|;|,"); if ( sqrt((a[5] - ((a[3] + 24) + 24))^2) >= 50 && a[6] == "+") print $1"\t"a[2]"\t"a[5]-50"\t"a[5]"\t"a[6]}' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header > ${cancer}_featureCounts_upstream_window
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,"=|:|;|,"); if ( sqrt((a[5] - ((a[3] + 24) + 24))^2) >= 50 && a[6] == "+") print $1"\t"a[2]"\t"a[5]"\t"a[5]+50"\t"a[6]}' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header > ${cancer}_featureCounts_downstream_window
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,"=|:|;|,"); if ( sqrt((a[5] - ((a[3] + 24) + 24))^2) < 50 && a[6] == "+") print $1"\t"a[2]"\t"a[5]-sqrt((a[5] - ((a[3] + 24) + 24))^2)"\t"a[5]"\t"a[6]}' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header >> ${cancer}_featureCounts_upstream_window
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,"=|:|;|,"); if ( sqrt((a[5] - ((a[3] + 24) + 24))^2) < 50 && a[6] == "+") print $1"\t"a[2]"\t"a[5]"\t"a[5]+sqrt((a[5] - ((a[3] + 24) + 24))^2)"\t"a[6]}' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header >> ${cancer}_featureCounts_downstream_window
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,"=|:|;|,"); if ( sqrt((a[5] - ((a[3] - 24) - 24))^2) >= 50 && a[6] == "-") print $1"\t"a[2]"\t"a[5]-50"\t"a[5]"\t"a[6]}' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header >> ${cancer}_featureCounts_upstream_window
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,"=|:|;|,"); if ( sqrt((a[5] - ((a[3] - 24) - 24))^2) >= 50 && a[6] == "-") print $1"\t"a[2]"\t"a[5]"\t"a[5]+50"\t"a[6]}' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header >> ${cancer}_featureCounts_downstream_window
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,"=|:|;|,"); if ( sqrt((a[5] - ((a[3] - 24) - 24))^2) < 50 && a[6] == "-") print $1"\t"a[2]"\t"a[5]-sqrt((a[5] - ((a[3] - 24) - 24))^2)"\t"a[5]"\t"a[6]}' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header >> ${cancer}_featureCounts_upstream_window
awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,"=|:|;|,"); if ( sqrt((a[5] - ((a[3] - 24) - 24))^2) < 50 && a[6] == "-") print $1"\t"a[2]"\t"a[5]"\t"a[5]+sqrt((a[5] - ((a[3] - 24) - 24))^2)"\t"a[6]}' peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq_header >> ${cancer}_featureCounts_downstream_window
grep -v "GeneID" ${cancer}_featureCounts_upstream_window | awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {print $0}' > ${cancer}_featureCounts_upstream_window1
grep -v "GeneID" ${cancer}_featureCounts_downstream_window | awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {print $0}' > ${cancer}_featureCounts_downstream_window1

for file in $(cat ${files}); do
        id=$(echo "${file}" | cut -d "." -f1)
	bam_file="${bam_path}/"*"bam"
        N=24
        if (( i % N == 0 )); then
                wait
        fi
        ((i++))
	(
	featureCounts -F SAF -p -s 0 -g gene_id -t exon -a ${cancer}_featureCounts_upstream_window1 -o ${id}_featureCounts_upstream_window_counts.txt ${bam_file} -T 6 -f -M -O --ignoreDup
	featureCounts -F SAF -p -s 0 -g gene_id -t exon -a ${cancer}_featureCounts_downstream_window1 -o ${id}_featureCounts_downstream_window_counts.txt ${bam_file} -T 6 -f -M -O --ignoreDup
	paste -d "\t" <(tail -n+3 ${id}_featureCounts_upstream_window_counts.txt | awk '{print $1}') <(tail -n+3 ${id}_featureCounts_upstream_window_counts.txt | awk '{print $NF}') <(tail -n+3 ${id}_featureCounts_downstream_window_counts.txt | awk '{print $NF}') | awk '{split($1,a,"=|:|;|,"); if (a[6] == "+" && $(NF-1)-$NF > 8 && $(NF-1) > 8 && $(NF-1)/($NF + $(NF-1)) > 2/3) print $1}' > ${id}_peptideSeqsFASTA_header_passed
	paste -d "\t" <(tail -n+3 ${id}_featureCounts_upstream_window_counts.txt | awk '{print $1}') <(tail -n+3 ${id}_featureCounts_upstream_window_counts.txt | awk '{print $NF}') <(tail -n+3 ${id}_featureCounts_downstream_window_counts.txt | awk '{print $NF}') | awk '{split($1,a,"=|:|;|,"); if (a[6] == "-" && $NF-$(NF-1) > 8 && $NF > 8 && $NF/($(NF-1) + $NF) > 2/3) print $1}' >> ${id}_peptideSeqsFASTA_header_passed
	awk '{split($1,a,":|;|,"); if (sqrt((a[4]-a[2])^2) > 48) print }' ${id}_peptideSeqsFASTA_header_passed > ${id}_peptideSeqsFASTA_header_passed1 ## iPASs will be deleted if the distance to intron boundary is no more than 24
	for entry in $(cat ${id}_peptideSeqsFASTA_header_passed1) ; do grep --no-group-separator "${entry}"$ peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq.fa -A 1 >> peptideSeqsFASTA_with_SigDrop_${id}.fa ; done
	) &
done

#for file in $(cat ${files}) ; do
#        id=$(echo "${file}" | cut -d "." -f1)
#        bedtools intersect -wa -wb -v -a <( awk '{split($1,a,/=|:|;|,/); print a[2]"\t"a[5]"\t"a[5]"\t"$1"\t0\t"a[6]}' ${id}_peptideSeqsFASTA_header_passed1 ) -b ${path}/IPA_events_from_GTEx_normal_for_filtering > ${id}_peptideSeqsFASTA_header_passed2
#        bedtools intersect -wa -wb -v -a ${id}_peptideSeqsFASTA_header_passed2 -b ${path}/IPA_events_from_tcga_normal_for_filtering > ${id}_peptideSeqsFASTA_header_passed3
#        bedtools intersect -wa -wb -v -a ${id}_peptideSeqsFASTA_header_passed3 -b ${path}/IPA_events_from_blueprint_normal_for_filtering > ${id}_peptideSeqsFASTA_header_passed4
#done

## intronic neo-epitope prediction
for file in $(cat ${files}); do
	id=$(echo "${file}" | cut -d "." -f1)
	bam_file="${bam_path}/${id}"*"bam"
	N=24
        if (( i % N == 0 )); then
                wait
        fi
        ((i++))
	(
	{ 
        read 
                while read -a line1; do 
            	        if [ ${line1[1]} == ${line1[3]} ]; then 
                	        hla1+="HLA-${line1[1]}," 
                        else 
       		        	if [ ${line1[1]} == "no" ] && [ ${line1[3]} != "no" ]; then
					hla1+="HLA-${line1[3]},"
				else
				        if [ ${line1[1]} != "no" ] && [ ${line1[3]} == "no" ]; then
						hla1+="HLA-${line1[1]},"
					else
						hla1+="HLA-${line1[1]},HLA-${line1[3]}," 
					fi
				fi
               	        fi 
                done  
	} < "${file_path}/${id}/${id}-ClassI-class.HLAgenotype4digits"
	
	for entry in $(cut -f4 ${id}_peptideSeqsFASTA_header_passed4) ; do grep --no-group-separator "${entry}"$ peptideSeqsFASTA_${cancer}_input_for_peptideseqs_uniq.fa -A 1 >> peptideSeqsFASTA_with_SigDrop_${id}.2.fa ; done
	awk '{ if (NR%2 != 0) line=$0; else {printf("%s\t%s\n", line, $0); line="";} } END {if (length(line)) print line;}' peptideSeqsFASTA_with_SigDrop_${id}.2.fa | awk -F'\t' 'NR>0{$0=$0"\t"NR-0} 1' | awk '{print ">s="$3"\t"$2}' | xargs -n1 > peptideSeqsFASTA_input_of_${id}
	awk '{ if (NR%2 != 0) line=$0; else {printf("%s\t%s\n", line, $0); line="";} } END {if (length(line)) print line;}' peptideSeqsFASTA_with_SigDrop_${id}.2.fa | awk -F'\t' 'NR>0{$0=$0"\t"NR-0} 1' | awk '{print "s_"$3"\t"substr($1,4,length($1))}' | sed 's/\tr/\tchr/g' > key_file_of_${id}
	"${netMHCpan}" -f peptideSeqsFASTA_input_of_${id} -inptype 0 -l ${window} -s -xls -xlsfile ${file_path}${id}_NETMHCpan_out.xls -a ${hla1} -BA > "/data/11000039/e0149673/scratch/Projects/TCGA_unsorted_bam/intron_polyadenylated_peptide_for_MS/logs/${id}"
	awk '{print $3"\t"$0}' ${id}_NETMHCpan_out.xls > tmp_${id}
	awk -v OFS='\t' 'NR==FNR {h[$1] = $2;next} {print h[$1],$0}' key_file_of_${id} tmp_${id} | awk '{if ($NF != 0) print}' | tail -n+3 | grep -v ',CD99,' > ${id}_output
	
##count the expression level of peptide sequences and remove the peptides with low coverage
        awk '{split($1,a,":|;|,"); if (a[5] == "+") printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\n", $0, a[2]+$3*3-3, a[2]+$3*3-3+2, a[2]+$3*3+length($4)*3-4-2, length($4)*3+$3*3+a[2]-4}' ${id}_output > ${id}_output_+
        awk '{split($1,a,":|;|,"); if (a[5] == "-") printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\n", $0, a[2]-$3*3-length($4)*3+4, a[2]-$3*3-length($4)*3+4+2, a[2]-$3*3+3-2, a[2]-$3*3+3}' ${id}_output > ${id}_output_-
        cat ${id}_output_- >> ${id}_output_+
        awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,":|;|,"); print $1"\t"a[1]"\t"$(NF-3)"\t"$(NF-2)"\t"a[5]}' ${id}_output_+ > ${id}_featureCounts_first_half_sequence
        awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,":|;|,"); print $1"\t"a[1]"\t"$(NF-1)"\t"$NF"\t"a[5]}' ${id}_output_+ > ${id}_featureCounts_second_half_sequence
        featureCounts -F SAF -p -s 0 -g gene_id -t exon -a ${id}_featureCounts_first_half_sequence -o ${id}_featureCounts_first_half_sequence_counts.txt ${bam_file} -T 6 -f -M -O
        featureCounts -F SAF -p -s 0 -g gene_id -t exon -a ${id}_featureCounts_second_half_sequence -o ${id}_featureCounts_second_half_sequence_counts.txt ${bam_file} -T 6 -f -M -O
        paste -d "\t" ${id}_output_+ <(tail -n+3 ${id}_featureCounts_first_half_sequence_counts.txt | awk '{print $NF}') <(tail -n+3 ${id}_featureCounts_second_half_sequence_counts.txt | awk '{print $NF}') | awk '{if ($NF >= 5 && $(NF-1) >= 5) print }'> ${id}_output
	awk -F'\t' -v OFS='\t' '{split($1,a,","); print a[1],a[2],a[4]"\t"a[3]"\t"$4}' ${id}_output | awk -F"\t" '!seen[$5]++' > ${id}.reliables
	
## filter with uniprot
	bash ${path}/construct_fasta.sh ${id}.reliables ${cancer} ${path}
        java -jar ${path}/PeptideMatchCMD_1.0.jar -a query -i ${path}/sprot_index_human/ -Q ${file_path}${id}.reliables.fa -e -o ${id}.reliables.fa.out
        for line in $(grep 'No match' ${id}.reliables.fa.out | cut -f1) ; do grep --no-group-separator "${line}"$ ${id}.reliables >> ${id}.reliables.filtered_by_uniprot ; done
	) &
done

#rm *.0*
