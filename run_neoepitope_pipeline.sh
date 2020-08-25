#!/bin/bash
#PBS -N BLCA

:<<USE
1, call cluster 
2, calculate IPA coverage 
3, predict IPA fasta files 
4. predict IPA neoantigens
USE

installDIR=`grep installDIR run.configure | sed s/.*\=//`
window=`grep peptide_length run.configure | sed s/.*\=//`
cancer=`grep dataset_name run.configure | sed s/.*\=//`
bam_path=`grep bam_dir run.configure | sed s/.*\=//`
OUTDIR=`grep output_dir run.configure | sed s/.*\=//`
soft_clip_path=`grep soft_clip_path run.configure | sed s/.*\=//`

mkdir $OUTDIR
cd $installDIR
cp findPeak.cluster.pl C.EncodeGencode_merged_intron_selected_excluded_exon_- B.EncodeGencode_merged_intron_selected_excluded_exon_+ callCluster.pl gencode_PASs_extended_by_100_filtered A.hg19_intron.excluded.exon.bed GeneratePeptide.py $OUTDIR

files="${installDIR}/sample_lists/TCGA_${cancer}.list"
generate_table="${installDIR}/generate_ipa_table.R"
rscript=`grep rscript run.configure | sed s/.*\=//`
python=`grep py2 run.configure | sed s/.*\=//`
netMHCpan=`grep netMHCpan run.configure | sed s/.*\=//`
featureCounts=`grep featureCounts run.configure | sed s/.*\=//`

cd $OUTDIR
mkdir "$OUTDIR/logs"

for file in $(cat ${files}); do

	BAM="${bam_path}/${file}/"*".bam"
	JUNCBED="${bam_path}/${file}/"*"SJ.out.tab"
	soft_clip_file="${soft_clip_path}/${cancer}/${file}.unsorted.bam_out_PAS_sites_from_softclipping_based_on_S.3_As.txt"
	STAR_final_output="${bam_path}/${file}/"*".Log.final.out"
	read_length=$(grep "Average input read length" $STAR_final_output | awk '{print $NF}')

### call cluster 
	echo call cluster for all soft clipped reads
	awk '$4 !~ /N/' ${soft_clip_file} >> $OUTDIR/${file}.0 ##remove splicing reads
	awk 'NR==FNR{a[$2]++; next} a[$2]>=2' <(sort -k3,3r -k1,1 -k2,2n $OUTDIR/${file}.0) <(sort -k3,3r -k1,1 -k2,2n $OUTDIR/${file}.0) > $OUTDIR/${file}.0.1  # print the PASs that supported by at least two soft-clipping reads
	sort -k3,3r -k1,1 -k2,2n $OUTDIR/${file}.0.1 | awk '{if ($3 == "+") print $1"\t"$2"\t"$2+1"\t"$4"\t0\t"$3}' > $OUTDIR/${file}.softclip.polya.bed
	sort -k3,3r -k1,1 -k2,2n $OUTDIR/${file}.0.1 | awk '{if ($3 == "-") print $1"\t"$2-1"\t"$2"\t"$4"\t0\t"$3}' >> $OUTDIR/${file}.softclip.polya.bed
 
	perl callCluster.pl $OUTDIR/${file}.softclip.polya.bed $OUTDIR/${file}.softclip.polya.cluster $OUTDIR/${file}.softclip.polya.cluster.length 

done

## find peak 
echo find peak for all samples
for file in $(cat ${files}); do awk '{if ($6 == "+") print $1"\t"$2"\t"$3"\t"$4":"FILENAME"\t"$5"\t"$6;else print $1"\t"$2"\t"$3"\t"$4":"FILENAME"\t"$5"\t"$6}' ${file}.softclip.polya.cluster | sed 's/\.softclip.polya.cluster//g' | awk '{ if ($1 == "1" || $1 == "2" || $1 == "3" || $1 == "4" ||$1 =="5"||$1 == "6" ||$1 =="7" || $1 =="8"||$1=="9"||$1=="10"||$1=="11"||$1=="12"||$1=="13"||$1=="14"||$1=="15"||$1=="16"||$1=="17"||$1=="18"||$1=="19"||$1=="20"||$1=="21"||$1=="22"||$1=="X"||$1=="Y") print}' >> $OUTDIR/${cancer}.cluster.bed ; done
sort -k1,1 -k2,2n $OUTDIR/${cancer}.cluster.bed > $OUTDIR/${cancer}.cluster.sorted.bed
mergeBed -d 24 -s -i $OUTDIR/${cancer}.cluster.sorted.bed -c 4,5,6 -o collapse,sum,distinct > $OUTDIR/${cancer}.cluster.merged.bed
for file in $(cat ${files}); do cat $OUTDIR/${file}.softclip.polya.bed >> $OUTDIR/${cancer}_all_SC_reads ; done
perl findPeak.cluster.pl $OUTDIR/${cancer}.cluster.merged.bed $OUTDIR/${cancer}_all_SC_reads $OUTDIR/${cancer}.cluster.merged_withPeak.bed
awk -F ',' '{if (NF >= 2 ) print }' $OUTDIR/${cancer}.cluster.merged_withPeak.bed | awk '{print "chr"$1"\t"$7"\t"$8"\t"$4"\t"$5"\t"$6}' > $OUTDIR/${cancer}.cluster.merged_withPeak.reliable.bed ##only keep the polyA clusters that supported by at least 2 samples for analyzing one cancer. But if there is only very few samples to be analyzed, then change it to one.

## extract only intron
bedtools intersect -wa -wb -a $OUTDIR/${cancer}.cluster.merged_withPeak.reliable.bed -b A.hg19_intron.excluded.exon.bed | awk '{if (($6 == "+" && $12 == "+") || ($6 == "-" && $12 == "-")) print}' |awk -F"\t" '!seen[$1, $2, $3, $6]++' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $OUTDIR/${cancer}.cluster.merged_withPeak.reliable.intron.bed

## HLA types
for file in $(cat ${files}); do
	id=$(echo "${file}" | cut -d "." -f1)
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
         } < "${soft_clip_path}/hla_types_newest/${id}-ClassI-class.HLAgenotype4digits"
	echo ${hla1} > $OUTDIR/${id}.hla
	hla1=""
done

for file in $(cat ${files}); do

	BAM="${bam_path}/${file}/"*".bam"
	JUNCBED="${bam_path}/${file}/"*"SJ.out.tab"
	soft_clip_file="${soft_clip_path}/${cancer}/${file}.unsorted.bam_out_PAS_sites_from_softclipping_based_on_S.3_As.txt"
	STAR_final_output="${bam_path}/${file}/"*".Log.final.out"
	read_length=$(grep "Average input read length" $STAR_final_output | awk '{print $NF}')
	grep "Uniquely mapped reads number" $STAR_final_output | awk '{print $NF}' > $OUTDIR/${file}.total.aligned.all.count

	id=$(echo "${file}" | cut -d "." -f1)
	hla1=$(cat ${id}.hla)
#	N=24
#        if (( i % N == 0 )); then
#                wait
#        fi
#        ((i++))
#        (

	grep -oP "(2[0-2]|1[0-9]|[1-9]|X|Y)\|\D+NA\D\d+\|\d+:(?=${file})" $OUTDIR/${cancer}.cluster.merged_withPeak.reliable.intron.bed |sed 's/://g' > $OUTDIR/${file}.intron.clusterID
	for line in $(cat ${file}.intron.clusterID); do grep "${line}" $OUTDIR/${cancer}.cluster.merged_withPeak.reliable.intron.bed | awk -v var=${line} '{split(var,a,/\|/); if ($6 == "+") print $1"\t"a[4]"\t"a[4]+1"\t"var"\t"a[5]"\t"$6;else print $1"\t"a[4]-1"\t"a[4]"\t"var"\t"a[5]"\t"$6}' >> $OUTDIR/${file}.intron.clusterID.output ;done
	grep -v "${file}" $OUTDIR/${cancer}.cluster.merged_withPeak.reliable.intron.bed | awk '{if ($6 == "+") print $1"\t"$2"\t"$3"\t"$1"|"$6"|NA|"$2"|"0"\t"0"\t"$6;else print $1"\t"$2"\t"$3"\t"$1"|"$6"|NA|"$3"|"0"\t"0"\t"$6}' >> $OUTDIR/${file}.intron.clusterID.output


## extract polyA flanking up- and down-stream 48nt
	awk '{if ($6 == "+") print $1"\t"$2-24"\t"$2+24"\t"$4"\t"$5"\t"$6;else print $1"\t"$3-24"\t"$3+24"\t"$4"\t"$5"\t"$6}' $OUTDIR/${file}.intron.clusterID.output > $OUTDIR/${file}.softclip.polya.win48

### calculate coverage 
	echo calculate coverage 
	echo -e "GeneID\tChr\tStart\tEnd\tStrand\t" > $OUTDIR/${file}.softclip.polya.ups100.SAF
	echo -e "GeneID\tChr\tStart\tEnd\tStrand\t" > $OUTDIR/${file}.softclip.polya.dns100.SAF

	awk '{OFS="\t"} {if($6 == "+") print $4,$1,$2 - 100,$2,"+";else print $4,$1,$3,$3 + 100,"-"}' $OUTDIR/${file}.softclip.polya.win48 >> $OUTDIR/${file}.softclip.polya.ups100.SAF
	awk '{OFS="\t"} {if($6 == "+") print $4,$1,$3,$3 + 100,"+";else print $4,$1,$2 - 100,$2,"-"}' $OUTDIR/${file}.softclip.polya.win48 >> $OUTDIR/${file}.softclip.polya.dns100.SAF

	$featureCounts -F SAF -p -s 0 -g gene_id -t exon -a $OUTDIR/${file}.softclip.polya.ups100.SAF -o $OUTDIR/${file}.softclip.polya.ups100.SAF.all.count $BAM -T 20 -f -M -O --ignoreDup
	$featureCounts -F SAF -p -s 0 -g gene_id -t exon -a $OUTDIR/${file}.softclip.polya.ups100.SAF -o $OUTDIR/${file}.softclip.polya.ups100.SAF.junc.count $BAM -T 20 -f -M -O --ignoreDup --splitOnly
	$featureCounts -F SAF -p -s 0 -g gene_id -t exon -a $OUTDIR/${file}.softclip.polya.dns100.SAF -o $OUTDIR/${file}.softclip.polya.dns100.SAF.all.count $BAM -T 20 -f -M -O --ignoreDup
	$featureCounts -F SAF -p -s 0 -g gene_id -t exon -a $OUTDIR/${file}.softclip.polya.dns100.SAF -o $OUTDIR/${file}.softclip.polya.dns100.SAF.junc.count $BAM -T 20 -f -M -O --ignoreDup --splitOnly

	awk '$4 == 1 || $4 == 3 || $4 == 5' $JUNCBED | awk '{OFS="\t"}{print $1,$2-1,$2,$1":"$2"-"$3":"$6,$7,"+"}' > $OUTDIR/${file}.junction.ss5.bed
	awk '$4 == 2 || $4 == 4 || $4 == 6' $JUNCBED | awk '{OFS="\t"}{print $1,$3-1,$3,$1":"$2"-"$3":"$6,$7,"-"}' >> $OUTDIR/${file}.junction.ss5.bed
	awk '$4 == 1 || $4 == 3 || $4 == 5' $JUNCBED | awk '{OFS="\t"}{print $1,$3,$3+1,$1":"$2"-"$3":"$6,$7,"+"}' > $OUTDIR/${file}.junction.ss3.bed
	awk '$4 == 2 || $4 == 4 || $4 == 6' $JUNCBED | awk '{OFS="\t"}{print $1,$2,$2+1,$1":"$2"-"$3":"$6,$7,"-"}' >> $OUTDIR/${file}.junction.ss3.bed

### overlap with junction 
	awk '{OFS="\t"} {if($6 == "+") print $1,$2,$2 + 24,$4,$5,"+";else print $1,$3 - 24,$3,$4,$5,"-"}' $OUTDIR/${file}.softclip.polya.win48 > $OUTDIR/${file}.softclip.polya.ups24 
	awk '{OFS="\t"} {if($6 == "+") print $1,$3 - 24,$3,$4,$5,"+";else print $1,$2,$2 + 24,$4,$5,"-"}' $OUTDIR/${file}.softclip.polya.win48 > $OUTDIR/${file}.softclip.polya.dns24 

	echo -e "Chr\tStart\tEnd\tName\tScore\tStrand" > $OUTDIR/${file}.softclip.polya.ups24.ss5.wo
echo -e "Chr\tStart\tEnd\tName\tScore\tStrand" > $OUTDIR/${file}.softclip.polya.dns24.ss3.wo
	intersectBed -a $OUTDIR/${file}.softclip.polya.ups24 -b $OUTDIR/${file}.junction.ss5.bed -s -u >> $OUTDIR/${file}.softclip.polya.ups24.ss5.wo
	intersectBed -a $OUTDIR/${file}.softclip.polya.dns24 -b $OUTDIR/${file}.junction.ss3.bed -s -u >> $OUTDIR/${file}.softclip.polya.dns24.ss3.wo

## filter with RPM and generate filtered tables
	$rscript $generate_table $OUTDIR ${file} ${read_length}

## filter with the putative IPA transcripts found from TCGA normal, GTEx, BLUEPRINT normal samples, and GENCODE
	bedtools intersect -wa -wb -v -a <(sort -k6,6r -k1,1 -k2,2n $OUTDIR/${file}.ipa.table) -b <(sort -k6,6r -k1,1 -k2,2n gencode_PASs ) > $OUTDIR/${file}.ipa.filtered1.table
	bedtools intersect -wa -wb -v -a $OUTDIR/${file}.ipa.filtered1.table -b $installDIR/IPA_events_from_GTEx_normal_for_filtering > $OUTDIR/${file}.ipa.filtered2.table
	bedtools intersect -wa -wb -v -a $OUTDIR/${file}.ipa.filtered2.table -b $installDIR/IPA_events_from_tcga_normal_for_filtering > $OUTDIR/${file}.ipa.filtered3.table
	bedtools intersect -wa -wb -v -a $OUTDIR/${file}.ipa.filtered3.table -b $installDIR/IPA_events_from_blueprint_normal_for_filtering > $OUTDIR/${file}.ipa.filtered.table

## generate IPA derived peptide sequences
	bedtools intersect -wa -wb -a $OUTDIR/${file}.ipa.filtered.table  -b  B.EncodeGencode_merged_intron_selected_excluded_exon_+ | awk '{if (($6 == "+" && $12 == "+") || ($6 == "-" && $12 == "-")) print}' |awk -F"\t" '!seen[$1, $2, $3, $6]++' > EncodeGencode_merged_intron_selected_+_overlapped_with_${file}
	bedtools intersect -wa -wb -a  $OUTDIR/${file}.ipa.filtered.table -b  C.EncodeGencode_merged_intron_selected_excluded_exon_- | awk '{if (($6 == "+" && $12 == "+") || ($6 == "-" && $12 == "-")) print}' |awk -F"\t" '!seen[$1, $2, $3, $6]++' > EncodeGencode_merged_intron_selected_-_overlapped_with_${file}
	awk -F'\t' -v OFS='\t' '{ if ($6 == "+") print $1,$2,$10,$13,$14,$15,$6,$5,$11,$4}' EncodeGencode_merged_intron_selected_+_overlapped_with_${file} > input_for_peptideseqs_of_1_${file}
	awk -F'\t' -v OFS='\t' '{ if ($6 == "-") print $1,$3,$10,$13,$14,$15,$6,$5,$11,$4}' EncodeGencode_merged_intron_selected_-_overlapped_with_${file} >> input_for_peptideseqs_of_1_${file}
	awk -F"\t" '!seen[$1, $2, $(NF-2), $(NF-1), $NF]++' input_for_peptideseqs_of_1_${file} > ${file}_input_for_peptideseqs_uniq

	${python} GeneratePeptide.py "${window}" "${OUTDIR}" ${file}_input_for_peptideseqs_uniq
	grep '>s=' peptideSeqsFASTA_${file}_input_for_peptideseqs_uniq.fa > peptideSeqsFASTA_${file}_input_for_peptideseqs_uniq_header
	awk '{split($1,a,":|;|,"); if (sqrt((a[4]-a[2])^2) >= 48) print }' peptideSeqsFASTA_${file}_input_for_peptideseqs_uniq_header > ${file}_peptideSeqsFASTA_header_passed1 ## iPASs will be deleted if the distance to intron boundary is no more than 24
	for entry in $(cat ${file}_peptideSeqsFASTA_header_passed1) ; do grep --no-group-separator "${entry}"$ peptideSeqsFASTA_${file}_input_for_peptideseqs_uniq.fa -A 1 >> peptideSeqsFASTA_with_SigDrop_${file}.fa ; done

## intronic neo-epitope prediction
	awk '{ if (NR%2 != 0) line=$0; else {printf("%s\t%s\n", line, $0); line="";} } END {if (length(line)) print line;}' $OUTDIR/peptideSeqsFASTA_with_SigDrop_${id}.fa | awk -F'\t' 'NR>0{$0=$0"\t"NR-0} 1' | awk '{print ">s="$3"\t"$2}' | xargs -n1 > $OUTDIR/peptideSeqsFASTA_input_of_${id}
        awk '{ if (NR%2 != 0) line=$0; else {printf("%s\t%s\n", line, $0); line="";} } END {if (length(line)) print line;}' $OUTDIR/peptideSeqsFASTA_with_SigDrop_${id}.fa | awk -F'\t' 'NR>0{$0=$0"\t"NR-0} 1' | awk '{print "s_"$3"\t"substr($1,4,length($1))}' | sed 's/\tr/\tchr/g' > $OUTDIR/key_file_of_${id}
        "${netMHCpan}" -f $OUTDIR/peptideSeqsFASTA_input_of_${id} -inptype 0 -l ${window} -s -xls -xlsfile $OUTDIR/${id}_NETMHCpan_out.xls -a ${hla1} -BA > "$OUTDIR/logs/${id}"
	awk '{print $3"\t"$0}' $OUTDIR/${id}_NETMHCpan_out.xls > $OUTDIR/tmp_${id}
        awk -v OFS='\t' 'NR==FNR {h[$1] = $2;next} {print h[$1],$0}' $OUTDIR/key_file_of_${id} $OUTDIR/tmp_${id} | awk '{if ($NF != 0) print}' | tail -n+3 | grep -v ',CD99,' > $OUTDIR/${id}_output

        ##count the expression level of peptide sequences
        awk '{split($1,a,":|;|,"); if (a[5] == "+") printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\n", $0, a[2]+$3*3-3, a[2]+$3*3-3+2, a[2]+$3*3+length($5)*3-4-2, length($5)*3+$3*3+a[2]-4}' $OUTDIR/${id}_output > $OUTDIR/${id}_output_+
        awk '{split($1,a,":|;|,"); if (a[5] == "-") printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\n", $0, a[2]-$3*3-length($5)*3+4, a[2]-$3*3-length($5)*3+4+2, a[2]-$3*3+3-2, a[2]-$3*3+3}' $OUTDIR/${id}_output > $OUTDIR/${id}_output_-
        cat $OUTDIR/${id}_output_- >>$OUTDIR/${id}_output_+
        awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,":|;|,"); print $1"\t"a[1]"\t"$(NF-3)"\t"$(NF-2)"\t"a[5]}' $OUTDIR/${id}_output_+ > $OUTDIR/${id}_featureCounts_first_half_sequence
        awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand"} {split($1,a,":|;|,"); print $1"\t"a[1]"\t"$(NF-1)"\t"$NF"\t"a[5]}' $OUTDIR/${id}_output_+ > $OUTDIR/${id}_featureCounts_second_half_sequence
        $featureCounts -F SAF -p -s 0 -g gene_id -t exon -a $OUTDIR/${id}_featureCounts_first_half_sequence -o $OUTDIR/${id}_featureCounts_first_half_sequence_counts.txt ${BAM} -T 6 -f -M -O
        $featureCounts -F SAF -p -s 0 -g gene_id -t exon -a $OUTDIR/${id}_featureCounts_second_half_sequence -o $OUTDIR/${id}_featureCounts_second_half_sequence_counts.txt ${BAM} -T 6 -f -M -O
        paste -d "\t" ${id}_output_+ <(tail -n+3 $OUTDIR/${id}_featureCounts_first_half_sequence_counts.txt | awk '{print $NF}') <(tail -n+3 $OUTDIR/${id}_featureCounts_second_half_sequence_counts.txt | awk '{print $NF}') | awk '{if ($NF >= 5 && $(NF-1) >= 5) print }'> $OUTDIR/${id}_output

        awk -F'\t' -v OFS='\t' '{split($1,a,","); print a[1],a[2],a[4]"\t"a[3]"\t"$4}' $OUTDIR/${id}_output | awk -F"\t" '!seen[$5]++' > $OUTDIR/${id}.reliables

	## filter with uniprot
        bash $installDIR/construct_fasta.sh $OUTDIR/${id}.reliables ${cancer} "${OUTDIR}"

        java -jar $installDIR/PeptideMatchCMD_1.0.jar -a query -i $installDIR/sprot_index_human/ -Q $OUTDIR/${id}.reliables.fa -e -o $OUTDIR/${id}.reliables.fa.out

        for line in $(grep 'No match' $OUTDIR/${id}.reliables.fa.out | cut -f1) ; do grep --no-group-separator "${line}"$ $OUTDIR/${id}.reliables >> $OUTDIR/${id}.reliables.filtered_by_uniprot ; done
#        ) &
done
