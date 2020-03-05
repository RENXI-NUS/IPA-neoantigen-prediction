# Intronic polyadenylated (IPA) neoantigen prediction

Synopsis:

This pipeline can predict the IPA-derived neoantigens based on large-scale RNA-seq datasets.


Usage:

In order to run, please:
1) Download and install bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html) and featureCounts (http://subread.sourceforge.net/).
2) Change the python (version 2) path in run_neoepitope_pipeline.sh file (line 12).
3) Download NetMHCPan-4.0 (http://www.cbs.dtu.dk/services/NetMHCpan/) and change paths in runNetMHCpan.py file (line 14).
4) Download twoBitToFa utility from UCSC genome browser (https://genome.ucsc.edu/goldenpath/help/twoBit.html) and change the path in GeneratePeptide.py file (line 173).
5) Download Seq2HLA (https://github.com/TRON-Bioinformatics/seq2HLA) and change the path in run_neoepitope_pipeline.sh file (line 106) after using it for HLA typing.
6) Download all the four reference files from https://drive.google.com/drive/folders/1otsrhoVvqL-XHElO4LI_pVd1Aij1CVrk?usp=sharing.
    A. splicing junction coordinates derived from patients of 12 TCGA cancers. It could be generated customizedly from STAR output *.sjdb.out.  
    B. BED fromat files containing all introns from GENCODE annotation.  
7) Input: bam files sorted by name within the same directory (Prefer output file from STAR)    
8) Start the pipeline with run.sh file.
