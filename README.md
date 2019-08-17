ChIA-PET Tool V3
=
Introduction
=
Chromatin Interaction Analysis with Paired-End Tag (ChIA-PET) sequencing is a technology to study genome-wide long-range chromatin interactions bound by protein factors. ChIA-PET Tool V3, a software package for automatic processing of ChIA-PET sequence data, including: 
1.	linker filtering
2.	mapping the paired-end reads to a reference genome
3.	purifying the mapped reads
4.	dividing the reads into different categories
5.	peak calling
6.	interaction calling
7.	visualizing the results

Install
=
Java is a popular platform-independent programming language and can be run on any machines with a Java Virtual Machine (JVM). BWA is used to map ChIA-PET sequencing reads to a reference genome. SAMtools is used to convert the alignment output from SAM format to BAM format. BEDTools is required to convert the files from BAM format to bed format. R environment and its packages are used to compute the p-values in peak calling and interaction calling and generate the graphs for visualization.
ChIA-PET Tool V3 requires the following softwares:
1.	JDK>=1.8(https://www.oracle.com/technetwork/java/javase/downloads/index.html)
2.	BWA(http://bio-bwa.sourceforge.net/)
3.	SAMtools(http://samtools.sourceforge.net/)
4.	BEDTools(https://bedtools.readthedocs.io/en/latest/)
5.	R(https://www.r-project.org/)
6.	R package grid(install.packages("grid"))
7.	R package xtable(install.packages("xtable"))
8.	R package RCircos(install.packages("RCircos"))

Download the ChIA-PET Tool V3 package from https://github.com/GuoliangLi-HZAU/ChIA-PET_Tool_V3. 
Unpack the package using the following command in your selected directory:  
`$ unzip ChIA-PET_Tool_V3.zip`

Test data sets
=
short-read test data set is part of RNAPII ChIA-PET data from human K562 cells: https://1drv.ms/u/s!AqzVTcWMvT40bhONKvBSCZDxqjA

long-read test data set is part of CTCF ChIA-PET data from human GM12878 cells: https://1drv.ms/u/s!AqzVTcWMvT40bzU5gKTuxxGRvz8

RNAPII ChIA-PET data from human K562 cells:<br>
  GSM832464 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM832464)
  
CTCF ChIA-PET data from human GM12878 cells:<br>
  GSM1872886 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1872886)

Usage
=
Before excuting the ChIA-PET Tool V3, you need to create genome index by BWA referring to http://bio-bwa.sourceforge.net/bwa.shtml and configure environment variables of bwa, samools and bamToBed (bedtools). After that, we can simply run it with one command line:

`$ java -jar ChIA-PET.jar [options]`

    The options are as follows:
    
    Necessary options:
    --mode: There are two modes for ChIA-PET Tool V3. 0 for short read, 1 for long read.
    --fastq1: path of read1 fastq file.
    --fastq2: path of read2 fastq file.
    --linker: linker file.
    --minimum_linker_alignment_score: Specifies the allowed minimum alignment score.
    --GENOME_INDEX: specifies the path of BWA index file.
    --GENOME_LENGTH: specifies the number of base pairs in the whole genome.
    --CHROM_SIZE_INFO: specifies the file that contains the length of each chromosome.
    --CYTOBAND_DATA: specifies the ideogram data used to plot intra-chromosomal peaks and interactions.
    --SPECIES: specifies the genome used to plot inter-chromosomal interactions, 1 for human, 2 for mouse and 3 for 
    others.
    
    Other options:
    --start_step: start with which step, 1: linker filtering; 2: mapping to genome; 3: removing redundancy; 4: 
    categorization of PETs; 5: peak calling; 6: interaction calling; 7: visualizing, default: 1"
    --output: specifies the directory to store the output data from ChIA-PET Tool V3, please specify the absolute path of output directory.
    --prefix: specifies the prefix of all the output files, default: out.
    --minimum_tag_length: Specifies the minimum tag length. Tag is the sequence after linker removal. This parameter 
    is better to be set  above 18bp. Default: 18.
    --maximum_tag_length: Specifies the maximum tag length. Default:1000 Specifies the maximum tag length. Default: 
    1000.
    --minSecondBestScoreDiff: Specifies the score difference between the best-aligned and the second-best aligned 
    linkers. Default: 3.
    --output_data_with_ambiguous_linker_info: Determines whether to print the linker-ambiguity PETs. 0: not print; 1: 
    print, Default: 1.
    --thread: the number of threads used in linker filtering and mapping to genome. Default: 1.
    --MAPPING_CUTOFF: The mapping threshold to remove the PETs with poor quality. Default: 30.
    --MERGE_DISTANCE: specifies the distance limit to merge the PETs with similar mapping locations. Default: 2.
    --SELF_LIGATION_CUFOFF: specifies the distance threshold between self-ligation PETs and intra- chromosomal 
    inter-ligation PETs. Default: 8000.
    --EXTENSION_LENGTH: specifies the extension length from the location of each tag, which is determined by the 
    median span of the self-ligation PETs. Default: 500.
    --MIN_COVERAGE_FOR_PEAK: specifies the minimum coverage to define peak regions. Default:5.
    --PEAK_MODE: There are two modes for peak calling. Number 1 is "peak region" mode, which takes all the 
    overlapping PET regions above the coverage threshold as peak regions, and number 2 is "peak summit" mode, which 
    takes the highest coverage of overlapping regions as peak regions. Default: 2.
    --MIN_DISTANCE_BETWEEN_PEAK: specifies the minimum distance between two peaks. If the distance of two peak 
    regions is below the threshold, only the one with higher coverage will be kept. Default: 500.
    --GENOME_COVERAGE_RATIO: specifies the estimated proportion of the genome covered by the reads. Default: 0.8.
    --PVALUE_CUTOFF_PEAK: specifies p-value to filter peaks that are not statistically significant. Default: 0.00001. 
    --INPUT_ANCHOR_FILE: a file which contains user-specified anchors for interaction calling. If you don't have this 
    file, please specify the value of this variable as "null" instead. Default: null.
    --PVALUE_CUTOFF_INTERACTION: specifies p-value to filter false positive interactions. Default:0.05.
    
Especially, the directories of data should be set properly to make sure that the programs could run smoothly. ChIA-PET Tool V3 will create a folder named by `OUTPUT_PREFIX` in the `OUTPUT_DIRECTORY`. The default value of `OUTPUT_DIRECTORY` is in the master folder "`ChIA-PET_Tool_V3/`" and `OUTPUT_PREFIX` is "out". Examples of `CHROM_SIZE_INFO` and `CYTOBAND_DATA` are both in the master folder "`ChIA-PET_Tool_V3/chromInfo/`". We recommend users to select one or create a new one in that folder.

The results will be visualized by a HTML file which is in the output folder "`OUTPUT_DIRECTORY/OUTPUT_PREFIX/OUTPUT_PREFIX.ChIA-PET_Report`".

Example running
=
short-read:<br>
java -jar ChIA-PET.jar --mode 0 --fastq1 test_short_1.fastq --fastq2 test_short_2.fastq --linker ChIA-PET_Tool_V3/linker/linker.txt --minimum_linker_alignment_score 8 --GENOME_INDEX hg19.fa --GENOME_LENGTH 3E9 --CHROM_SIZE_INFO ChIA-PET_Tool_V3/chromInfo/hg19.chromSize.txt --CYTOBAND_DATA ChIA-PET_Tool_V3/chromInfo/hg19_cytoBandIdeo.txt --SPECIES 1 --output test_short --prefix K562 --thread 4

long-read:<br>
java -jar ChIA-PET.jar --mode 1 --fastq1 test_long_1.fastq --fastq2 test_long_2.fastq --linker ChIA-PET_Tool_V3/linker/linker_long.txt --minimum_linker_alignment_score 14 --GENOME_INDEX hg19.fa --GENOME_LENGTH 3E9 --CHROM_SIZE_INFO ChIA-PET_Tool_V3/chromInfo/hg19.chromSize.txt --CYTOBAND_DATA ChIA-PET_Tool_V3/chromInfo/hg19_cytoBandIdeo.txt --SPECIES 1 --output test_long --prefix GM12878 --thread 4

Result file
=
* Example of peak file named `OUTPUT_PREFIX`.peak.FDRfiltered.txt

|chrom|summit start|summit end|peak coverage|p-value |p.adjust|
|-----|------------|----------|-------------|--------|--------|
|chr1 |840305      |840556    |21           |1.18e-08|2.57e-07|

chrom: chromosome name.

summit start: The start coordinate of peak summit.

summit end: The end coordinate of peak summit.

peak coverage: The highest coverage by tags in a peak.

p-value: This value represents the statistical significance of a peak, which is calculated by Poisson distribution.

p.adjust: P.adjust means p-value adjusted with Benjamini-Hockberg method for multiple hypothesis testing.

* Example of interaction file named `OUTPUT_PREFIX`.cluster.FDRfiltered.txt

|chrom1|start1|end1  |chrom2|start1|end2  |ipet counts|type|distance|tag count within anchor 1|tag count within anchor 2|p-value|p.adjust|-log10(p-value)|-log10(p.adjust)|
|------|------|------|------|------|------|-----------|----|--------|-------------------------|-------------------------|--------------------|--------|---------------|----------------|
|chr7  |562341|563231|chr7  |574231|575143|2          |1   |11901   |3                        |9                        |7.11e-11            |1.99E-10|10.15          |9.7             |

chrom1: The name of the chromosome on which the cluster anchor 1 exists.

start1: The start coordinate of cluster anchor 1.

end1: The end coordinate of cluster anchor 1.

chrom2: The name of the chromosome on which the cluster anchor 2 exists.

start2: The start coordinate of cluster anchor 2.

end2: The end coordinate of cluster anchor 2.

ipet count: Number of PETs between cluster anchor 1 and cluster anchor 2.

type: Interactions type. 1 represents intra-chromosomal interaction, and 0 represents inter-chromosomal interaction.

distance: Distance between anchors of an intra-chromosomal interaction cluster. If two anchors are located on different chromosomes, the value is set to 2,147,483,647.

tag count within anchor 1: Number of tags that fall in the cluster anchor 1.

tag count within anchor 2: Number of tags that fall in the cluster anchor 2.

p-value: This value represents the statistical significance of a chromatin interaction, which is calculated by hyper-geometric distribution.

p.adjust: P.adjust means p-value adjusted with Benjamini-Hockberg method for multiple hypothesis testing.

-log10(p-value): The negative logarithm of p-value.

-log10(p.adjust): The negative logarithm of adjusted p-value.

References
=
1.	Li G, Cai L, Chang H, Hong P, Zhou Q, Kulakova EV, Kolchanov NA, Ruan Y. Chromatin Interaction Analysis with Paired-End Tag (ChIA-PET) sequencing technology and application. BMC Genomics, 2014, 15, S11
2.	Li G, Chen Y, Snyder MP, Zhang MQ. ChIA-PET2: a versatile and flexible pipeline for ChIA-PET data analysis. Nucleic Acids Research, 2016: 1-10
3.	Li G, Fullwood MJ, Xu H, Mulawadi FH, Velkov S, Vega V, Ariyaratne PN, Mohamed YB, Ooi H, Tennakoon C, Wei C, Ruan Y, Sung W. Software ChIA-PET tool for comprehensive chromatin interaction analysis with paired-end tag sequencing. Genome Biology 2010, 11, R22
4.	Li X, Luo OJ, Wang P, Zheng M, Wang D, Piecuch E, Zhu JJ, Tian SZ, Tang Z, Li G, Ruan Y. Long-read ChIA-PET for base-pair-resolution mapping of haplotype-specific chromatin interactions. Protocol 2017, 12: 899-915
5.	Paulsen J, Rødland EA, Holden L, Holden M, Hovig E. A statistical model of ChIA-PET data for accurate detection of chromatin 3D interactions. Nucleic Acids Research, 2014, 42, e143
6.	Phanstiel DH, Boyle AP, Heidari N, Snyder MP. Mango: a bias-correcting ChIA-PET analysis pipeline. Genome analysis, 2015, 31: 3092-3098

Contact
=
guoliang.li@mail.hzau.edu.cn
