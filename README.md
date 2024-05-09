# ClusteringAnalysis_ermB_2024
 
Machine learning-based classification reveals distinct clusters of non-coding genomic allelic variations associated with Erm-mediated antibiotic resistance

Yongjun Tan, Alexandre Le Scornet, Mee-Ngan Frances Yap, Dapeng Zhang

This repository is to store the data and source code for the manuscript entitled "Machine learning-based classification reveals distinct clusters of non-coding genomic allelic variations associated with Erm-mediated antibiotic resistance".

Please cite the paper:
Tan Y, Scornet AL, Yap M-NF, Zhang D. Machine learning classification reveals distinct clusters of non-coding genomic allelic variations associated with antibiotic resistance. (In preparation)

The repository consists of two main parts.

Part One: Allele frequency calculation
1. ermB.align.list
   This file is the input file of ermB_frequency.py and contains a multiple sequence alignment (MSA) with a total of 21,525 sequences used in this study. The MSA was generated using BLASTN (output format 1).
2. ermB.align.table
   This file is the output plaintext table of ermB_frequency.py. User could open it using Notepad software in Windows. It contains several columns: the first column represents nucleotides for each position, the second column indicates the position, the third column shows the percentage of allele frequency based on the number of aligned sequences, the fourth column displays the number of aligned sequences at each position, the fifth column shows the percentage of allele frequency based on the total number of sequences, and the sixth column represents the total number of sequences. Additional columns from the seventh onward categorize and count variations at each position.
4. ermB_frequency.py
   This file is the Python3 scripts to calculate allele frequency. Input file is ermB.align.list, and output file is ermB.align.table.

Part Three: Entropy calculation
1. ermB.align_1-211.list
   This file is the input file of ermB_clustering.py. This file was manually extracted from ermB.align.list. Because, we specifically focus on investigating the 1-211 region of upstream of the ermB.
2. kmeans_211.csv
   This file is the output file of ermB_clustering.py. We included cluster assignment results for 21,525 sequences from a series of trials with cluster numbers ranging from 20 to 100. The first row denotes each trial with a different cluster number, spanning from 20 to 100. The subsequent rows display the cluster assignments for each sequence, maintaining the same order as in the input file.
5. ermB_clustering.py
   This file is the Python3 scripts to conduct clustering analysis for the 1-211 region of upstream of the ermB. Input file is ermB.align_1-211.list, output file is kmeans_211.csv. In this analysis, we explored the optimal number of clusters for the K-means algorithm by conducting a series of trials with cluster numbers ranging from 20 to 100. All results are saved in the output file for subsequent analysis.


Comments, suggestions and bug reports to: yongjun.tan@slu.edu

