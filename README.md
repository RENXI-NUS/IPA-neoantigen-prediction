# Intronic polyadenylated (IPA) neoantigen prediction

**Synopsis:**

This pipeline can predict the IPA-derived neoantigens based on large-scale RNA-seq datasets.



**Usage:**

In order to run, please:
1) Download and install bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html), featureCounts (http://subread.sourceforge.net/), NetMHCPan-4.0 (http://www.cbs.dtu.dk/services/NetMHCpan/) and Seq2HLA (https://github.com/TRON-Bioinformatics/seq2HLA).
2) Download twoBitToFa and faToTwoBit utilities from UCSC genome browser (https://genome.ucsc.edu/goldenpath/help/twoBit.html) or install directly from conda.
3) Download all the four reference files from https://drive.google.com/drive/folders/1otsrhoVvqL-XHElO4LI_pVd1Aij1CVrk?usp=sharing, and put them into the IPA-neoantigen-prediction folder.

    A: Splicing junction coordinates derived from patients of 12 TCGA cancers. This file could be customizedly generated and substituted by processing the STAR output files ending *.sj.out.tab    
    B: BED fromat file containing all introns from GENCODE annotation    
    C: File of intron boundaries with frame information from positive strand    
    D: File of intron boundaries with frame information from negative strand
4) Download reference files "_c_Lucene41_0.pos" and "_c_Lucene41_0.doc" from https://drive.google.com/drive/folders/1otsrhoVvqL-XHElO4LI_pVd1Aij1CVrk?usp=sharing, and put them into the IPA-neoantigen-prediction/sprot_index_human folder
5) Configure the run.configure file according to the instruction of each input.
6) Ensure all the scripts are executable (chmod +x *) and start the pipeline with run.sh file.



**System Requirements:**

1) OS Requirements:
    Linux Ubuntu 16.04
2) Python (tested on v2.7.17) dependencies (install the latest versions available):
    numpy
    subprocess
    Bio.Seq
3) R (tested on v4.1.2) dependencies (install the latest versions available): 
    dplyr
    caroline
4) Perl (tested on v5.26.2)
5) Software tested version (install tested version: conda install -c bioconda bedtools=2.23.0):

    bedtools v2.23.0
    
    samtools v0.1.19
    
    STAR v2.7.10a
    
    NetMHCpan v4.0
    
    Seq2HLA v2.2
       
    featureCounts v2.0.1
