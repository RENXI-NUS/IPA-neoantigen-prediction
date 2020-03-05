#!/bin/bash

##provide full path of current directory (STRING without space; use _ for space):
curr_dir=/data/11000039/e0149673/scratch/Projects/TCGA_unsorted_bam/scripts/test/IPA-neoantigen-prediction

##provide peptide length (The length of peptide to be used for binding affinity prediction, and must be from 8 to 14 amino acids)
peptide_length=9

##dataset name (the dataset name of this prediction):
dataset_name=test

##full path for the bam files and raw fastq files (The corresponding full path for the unsorted bam files (i.e. sort by name) and raw fastq files of the analyzed dataset. STRING without space; use _ for space):
bam_dir=/data/11000039/e0149673/scratch/Projects/TCGA_unsorted_bam/scripts/test/IPA-neoantigen-prediction

##outputfile path (STRING without space):
output_dir=/data/11000039/e0149673/scratch/Projects/TCGA_unsorted_bam/scripts/test/IPA-neoantigen-prediction
