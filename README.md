# IPA-neoantigen-prediction

This pipeline can predict the IPA derived neoantigens based on large-scale paired-end RNA-seq data.

In order to run, please:
1) Download and install bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html) and featureCounts (http://subread.sourceforge.net/).
2) Change the python path in run_neoepitope_pipeline.sh file (line 12).
3) Download NetMHCPan-4.0 (http://www.cbs.dtu.dk/services/NetMHCpan/) and change paths in runNetMHCpan.py file (line 14).
4) Download twoBitToFa utility from UCSC genome browser (https://genome.ucsc.edu/goldenpath/help/twoBit.html) and change the path in PASRA2PeptideSeqs.py file (line 173).
5) Download Seq2HLA (https://github.com/TRON-Bioinformatics/seq2HLA) and change the path in run_neoepitope_pipeline.sh file (line 106) after using it for HLA typing.
