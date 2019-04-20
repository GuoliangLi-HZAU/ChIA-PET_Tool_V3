ChIA-PET Tool
===
Introduction
===
Chromatin Interaction Analysis with Paired-End Tag (ChIA-PET) sequencing is a technology to study genome-wide long-range chromatin interactions bound by protein factors. ChIA-PET Tool, a software package for automatic processing of ChIA-PET sequence data, including: 
1.	linker filtering
2.	mapping the paired-end reads to a reference genome
3.	purifying the mapped reads
4.	dividing the reads into different categories
5.	peak calling
6.	interaction calling
7.	visualizing the results

Install
===
Java is a popular platform-independent programming language and can be run on any machines with a Java Virtual Machine (JVM). BWA is used to map ChIA-PET sequencing reads to a reference genome. SAMtools is used to convert the alignment output from SAM format to BAM format. BEDTools is required to convert the files from BAM format to bed format. R environment and its packages are used to compute the p-values in peak calling and interaction calling and generate the graphs for visualization.
ChIA-PET Tool requires the following softwares:
1.	JDK>=1.8(https://www.oracle.com/technetwork/java/javase/downloads/index.html)
2.	BWA(http://bio-bwa.sourceforge.net/)
3.	SAMtools(http://samtools.sourceforge.net/)
4.	BEDTools(https://bedtools.readthedocs.io/en/latest/)
5.	R(https://www.r-project.org/)
6.  R package grid(install.packages("grid"))
7.	R package xtable(install.packages("xtable"))
8.	R package RCircos(install.packages("RCircos"))

Download the ChIA-PET Tool package from ***  
Unpack the package using the following command in your selected directory:  
`$ unzip ChIA_PET_Tool.zip`

Usage
===
Before excuting the ChIA-PET Tool, you need to create genome index by BWA referring to http://bio-bwa.sourceforge.net/bwa.shtml.
Edit the variables in the parameter file to set up the information about the data and parameters required for ChIA-PET Tool. Especially, the directories of data should be set properly to make sure that the programs could run smoothly. ChIA-PET Tool will create a folder named by `OUTPUT_PREFIX` in the `OUTPUT_DIRECTORY`. The default value of `OUTPUT_DIRECTORY` is in the master folder “`ChIA_PET_Tool/`”, and `OUTPUT_PREFIX` is “output”. We have to make sure the folder “`OUTPUT_DIRECTORY/OUTPUT_PREFIX`” is empty to avoid overwriting old files. The linker_file is in the master folder “`ChIA_PET_Tool /linker/`”. `CHROM_SIZE_INFO` and `CYTOBAND_DATA` are both in the master folder “`ChIA_PET_Tool /chromInfo/`”. You need to select one or create a new one in that folder.

    The parameters are as follows:
    OUTPUT_DIRECTORY: specifies the directory to store the output data from ChIA-PET Tool.
    OUTPUT_PREFIX: specifies the prefix of all the output files.
    MODE: There are two modes for ChIA-PET Tool. 0 for short read, 1 for long read. Default:0
    Fastq_file_1: The path and file name of read1 fastq file in pair-end sequencing.
    Fastq_file_2: The path and file name of read2 fastq file in pair-end sequencing.
    linker_file: A file specifies linker information, including linker sequences, the start positions and length of barcodes.
    minimum_linker_alignment_score: Specifies the allowed minimum alignment score. Default:8
    minimum_tag_length: Specifies the minimum tag length. Tag is the sequence after linker removal. This parameter is better to be set above 18bp. Default:20
    maximum_tag_length: Specifies the maximum tag length. Default:1000 Specifies the maximum tag length. Default:1000
    minSecondBestScoreDiff: Specifies the score difference between the best-aligned and the second-best aligned linkers. Default:3
    output_data_with_ambiguous_linker_info: Determines whether to print the linker-ambiguity PETs. Value 1 means to print these reads to specific files, while other positive integers mean not to print the PETs with ambiguous linkers. Default:1
    NTHREADS: the number of threads used in linker filtering and mapping. Default:4
    GENOME_INDEX: specifies the path of BWA index file.
    MAPPING_CUTOFF: The mapping threshold to remove the PETs with poor quality. Default:30
    MERGE_DISTANCE: specifies the distance limit to merge the PETs with similar mapping locations. Default:2
    SELF_LIGATION_CUFOFF: specifies the distance threshold between self-ligation PETs and intra- chromosomal inter-ligation PETs. Default:8000
    EXTENSION_LENGTH: specifies the extension length from the location of each tag, which is determined by the median span of the self-ligation PETs. Default:500
    MIN_COVERAGE_FOR_PEAK: specifies the minimum coverage to define peak regions. Default:5
    PEAK_MODE: There are two modes for peak calling. Number 1 is “peak region” mode, which takes all the overlapping PET regions above the coverage threshold as peak regions, and number 2 is “peak summit” mode, which takes the highest coverage of overlapping regions as peak regions. Default:2
    MIN_DISTANCE_BETWEEN_PEAK: specifies the minimum distance between two peaks. If the distance of two peak regions is below the threshold, only the one with higher coverage will be kept. Default:500
    GENOME_LENGTH: specifies the number of base pairs in the whole genome.
    GENOME_COVERAGE_RATIO: specifies the estimated proportion of the genome covered by the reads. Default:0.8
    PVALUE_CUTOFF_PEAK: specifies p-value to filter peaks that are not statistically significant. Default:0.00001
    CHROM_SIZE_INFO: specifies the file that contains the length of each chromosome.
    INPUT_ANCHOR_FILE: a file which contains user-specified anchors for interaction calling. If you don't have this file, please specify the value of this variable as “null” instead. Default:null
    PVALUE_CUTOFF_INTERACTION: specifies p-value to filter false positive interactions. Default:0.05
    CYTOBAND_DATA: specifies the ideogram data used to plot intra-chromosomal peaks and interactions.
    SPECIES`: specifies the genome used to plot inter-chromosomal interactions, 1 for human, 2 for mouse and 3 for others. Default:1
    
After setting all the required tools, data and parameters, you can simply run it with one command line:

`$ sh run.sh`

The results will be visualized by a HTML file which is in the output folder “`OUTPUT_DIRECTORY/OUTPUT_PREFIX/OUTPUT_PREFIX.ChIA-PET_Tool_Report`”.

References
===
1.	Li G, Cai L, Chang H, Hong P, Zhou Q, Kulakova EV, Kolchanov NA, Ruan Y. Chromatin Interaction Analysis with Paired-End Tag (ChIA-PET) sequencing technology and application. BMC Genomics, 2014, 15, S11
2.	Li G, Chen Y, Snyder MP, Zhang MQ. ChIA-PET2: a versatile and flexible pipeline for ChIA-PET data analysis. Nucleic Acids Research, 2016: 1-10
3.	Li G, Fullwood MJ, Xu H, Mulawadi FH, Velkov S, Vega V, Ariyaratne PN, Mohamed YB, Ooi H, Tennakoon C, Wei C, Ruan Y, Sung W. Software ChIA-PET tool for comprehensive chromatin interaction analysis with paired-end tag sequencing. Genome Biology 2010, 11, R22
4.	Li X, Luo OJ, Wang P, Zheng M, Wang D, Piecuch E, Zhu JJ, Tian SZ, Tang Z, Li G, Ruan Y. Long-read ChIA-PET for base-pair-resolution mapping of haplotype-specific chromatin interactions. Protocol 2017, 12: 899-915
5.	Paulsen J, Rødland EA, Holden L, Holden M, Hovig E. A statistical model of ChIA-PET data for accurate detection of chromatin 3D interactions. Nucleic Acids Research, 2014, 42, e143
6.	Phanstiel DH, Boyle AP, Heidari N, Snyder MP. Mango: a bias-correcting ChIA-PET analysis pipeline. Genome analysis, 2015, 31: 3092-3098

Contact
===

