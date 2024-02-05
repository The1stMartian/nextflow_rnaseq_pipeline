# RNA-Seq Pipeline
<i>An RNA-Seq pipeline for eukaryotic systems written in Python</i>
C. Breuer 2024

# Overview:
This pipeline is designed for the analysis of RNA-Seq data from eukaryotic systems. The offers quality assessment with FastQC, and calls pipeline RNA STAR, an efficient splicing-aware aligner via the linux (or WSL) command line. Additional scripts are provided for running differential expression analysis with DESeq2 in R, and for annotating those results with gene names using a GTF file. 

# Prerequisites
1) Linux or WSL. It is expected that this pipeline will be executed in a linux environment.
2) Python3 and Java (for FastQC)
- apt-install python3
3) The RNA-STAR aligner should be installed. This can be done using the following commands in Ubuntu:
- apt-get update
- apt-install RNA-STAR
- For more information, see the  [STAR manual].(https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).
4) An RNA-STAR reference genome (see the STAR manual for commands). This requires you to have a fasta formatted genome file and (optionally) an annotations file (e.g. a .gtf file). <b>It is highly beneficial to include the GTF annotation file as this will allow STAR to tally read counts while doing the alignmenment.</b> Place the output STAR genome files in a folder named "genomefiles".
5) For rDNA quality control mapping, an rRNA (or tRNA, etc.) fasta sequence and Bowtie2 are required. Bowtie2 can be installed with:
- apt-get update
- apt-install bowtie2
- Follow the instructions for installing a bowtie genome. e.g. bowtie2-build <path/to/genome.fa> <genomeName>. Move the output genome files to the "genomefiles" folder.


# Setup
1) Create a working folder for your analysis
2) Within the working folder, you should have your fastq files in a folder called "fastq". Fastq file name should be formatted like "sample1_R1.fastq" and "sample1_R2.fastq" for a typical paired-end sample. Files can also be gzipped and named like "Sample1_R1.fq.gz" or "Sample1_R1.fastq.gz". If additional formats are required, please edit the rnaSeq.py script, or change the filenames with the provided change_filenames.py script (no arguments needed).
3) Ensure that STAR can be called in the linux/WSL command line using command "RNA-STAR".
4) All STAR and Bowtie2 genome files should be placed in a subdirectory called "genomefiles". For example, for STAR genome named "dm6", the files should be under yourworkingdirectory/genomefiles/dm6. For bowtie2 genome "bs168", the files should be under yourworkingdirectory/genomefiles. Note that these filepaths can be changed if you edit the rnaSeq.py script customization section.
5) Fasta formatted genome files and genome annotation files (e.g. .gtf) should be in ./genomefiles. 
6) Edit the User Customization section of the rnaSeq.py script in any text editor. Specifically, change line 265 as needed to indicate the name of the genome you're mapping to. i.e. if you named the STAR genome "dm6" for Drosophila melanogaster build number 6, then that line should read g='dm6'.

# Running the pipeline
1) The pipeline can be run by executing the rnaSeq.py script with the command: python rnaSeq.py. 
2) Follow the prompts to activate FastQC, alignment by STAR, and to use FeatureCounts. 
3) Notes:
- Default settings will result in STAR producing read counts as part of the alignment process. This will only occur if you input a gene annotation file when building your STAR genome.
- STAR will output aligned files in .bam format. 