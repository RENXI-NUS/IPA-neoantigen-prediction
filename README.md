# Intronic polyadenylated (IPA) neoantigen prediction

Synopsis:

This pipeline can predict the IPA-derived neoantigens based on large-scale RNA-seq datasets.


Usage:

In order to run, please:
1) Download and install bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html), featureCounts (http://subread.sourceforge.net/), NetMHCPan-4.0 (http://www.cbs.dtu.dk/services/NetMHCpan/) and Seq2HLA (https://github.com/TRON-Bioinformatics/seq2HLA).
2) Download twoBitToFa utility from UCSC genome browser (https://genome.ucsc.edu/goldenpath/help/twoBit.html) and change the path in GeneratePeptide.py file (line 173).
3) Download all the four reference files from https://drive.google.com/drive/folders/1otsrhoVvqL-XHElO4LI_pVd1Aij1CVrk?usp=sharing, and put them into the IPA-neoantigen-prediction folder.

    A: Splicing junction coordinates derived from patients of 12 TCGA cancers. This file could be customizedly generated and substituted by processing the STAR output files ending *.sj.out.tab    
    B: BED fromat file containing all introns from GENCODE annotation    
    C: File of intron boundaries with frame information from positive strand    
    D: File of intron boundaries with frame information from negative strand
4) Download reference files "_c_Lucene41_0.pos" and "_c_Lucene41_0.doc" from https://drive.google.com/drive/folders/1otsrhoVvqL-XHElO4LI_pVd1Aij1CVrk?usp=sharing, and put them into the IPA-neoantigen-prediction/sprot_index_human folder
5) Input: 
    -  BAM files (sorted by name) are copyed within the same directory (output file from STAR is preferred and the name of the file should change to *.unsorted.bam)
    -  FASTQ files (for HLA typing) are copyed within the same directory (the naming of BAM file and corresponding FASTQ file should be consistent)
6) Configure the run.config file according to the instruction of each input.
7) Start the pipeline with run.sh file.
