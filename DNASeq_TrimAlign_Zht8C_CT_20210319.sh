#!/bin/bash

# folder with all folders of fastq files
top_level_folder="/data/guang/Experiment_20210312/raw_data/Zht8C_CT"


##******************* Parameter Definition Finished Here *********************##

##*************** The actual alignment script starts here. ********************##
################################################################################

## Global variables used by all script pieces.
## The global variables will be re-defined in each part. Although this re-definition
## is not necessary for execution of this merged script, re-definition of these
## global variables will make each part still a complete script and can be copied
## out to run independently.
## Re-define of the global variables to make  this part independent.
## top_level_folder="/data/guang/Experiment_20210312/raw_data/Zht8C_CT"
## Initialize the top_level folder. Must use full path!!
## This folder should contain the sample folders with single-end fatq.gz files.
## This folder will be used by all script pieces below.
cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"
## Get names of sample folders:
## The command on the right site will first use 'ls -l' to check all files and folders
## then 'grep "^d" ' will select the ones will a 'd' property at the beginning of 
## 'ls -l' command which are folders
## finally 'awk '{print $NF}' will print out the last column of 'ls -l' which are 
## the actual folder names. 
## NOTE: There should be no space in the folder names. If space exists in folder names
##       only the last word of folder name will be print out. This is not what we want!
## sample_folder_names will be used by all script pieces below.

################## Part 1 Quality Control and fastq data trimming.#####################
############# Part 1.1 Check fastq data quality using FastQC ##########################

for sample_folder in $sample_folder_names
do
	cd $top_level_folder/$sample_folder
	# Get into the sample_folder with fastq.gz file(s)
	mkdir -p $top_level_folder/$sample_folder/$sample_folder\_fastqc_results
	# Make a new folder in sample_folder to store FastQC result; Full path used here.
	fastqc -t 8 *.fq.gz -O ./$sample_folder\_fastqc_results/
	# -t 8: use 8 threads
	# relative path is used. The full path is 
	# $top_level_folder/$sample_folder/$sample_folder\_fastqc_results
	# Two files generate after calling fastqc: 
	#                (fastq.gz filename)_fastqc.html and (fastq.gz filename)_fastqc.zip
done


########################################################################################
################# Trimming with Trim_Galore
cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"

for sample_folder in $sample_folder_names
do
	cd $top_level_folder/$sample_folder
	# Get into the sample_folder with fastq.gz file
	mkdir -p $top_level_folder/$sample_folder/$sample_folder\_trim_galore_trimmed
	# Make a new folder in sample_folder to store trimmed result; Full path used here.
	fastq_gz_files="$(ls *.fq.gz)"
	# Vazyme use Nextera Transposase Adapters
	trim_galore --paired --retain_unpaired --phred33 --nextera \
	--length 20 \
	--output_dir ./$sample_folder\_trim_galore_trimmed \
	$fastq_gz_files
	
	# remove old fastq files to save space
	# rm *.gz
done


######### Parameters related to bwa-aln ###############
# bwa genome index is generated in TruSeq_totalRNA-Seq.sh script
genome_bwa_index_prefix="/data/guang/genome_index/mouse_bwa_index/GRCm38.101.dna.primary_assembly_plus_LambdaDNA.bwaIndex"
####################################################

top_level_folder="/data/guang/Experiment_20210312/raw_data/Zht8C_CT"

cd $top_level_folder
sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"

for sample_folder in $sample_folder_names
do
	cd $top_level_folder/$sample_folder
	#get into sample_folder
	mkdir -p $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results
	
	trimmed_fastq_gz_files="$(ls -d $PWD/$sample_folder\_*\_trimmed/*.fq.gz)"
	#$PWD is current path
	#This command get the full path of trimmed fast.gz files
	
	num_trimmed_fastq="$(echo "$trimmed_fastq_gz_files" | wc -l)"
	if [ $num_trimmed_fastq -ge 2 ] 
	# Check numer of trimmed fastq.gz files. If great than or equals to 2, they must come from 
	# paired-end read files.
	# Treat as trimmed fastq files from paired-end fastq files
	then
	paired_trimmed_read1="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_val_1.fq.gz)"
	paired_trimmed_read1_basename="$(basename $paired_trimmed_read1)"
	paired_trimmed_read2="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_val_2.fq.gz)"
	paired_trimmed_read2_basename="$(basename $paired_trimmed_read2)"
	unpaired_trimmed_reads="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_unpaired*.fq.gz)"
	# use unpaired_reads to store the filenames of the 2 unpaired read files trimmed out by trimmomatic
	
	# process trim galore trimmed read1 file
	bwa aln -t 8 $genome_bwa_index_prefix \
				 $paired_trimmed_read1 \
				 > $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results/$paired_trimmed_read1_basename.bwa.sai
	# process trim galore trimmed read1 file
	bwa aln -t 8 $genome_bwa_index_prefix \
				 $paired_trimmed_read2 \
				 > $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results/$paired_trimmed_read2_basename.bwa.sai
	# Get SAM file from paired-end trimmed fastq files
	bwa sampe $genome_bwa_index_prefix \
			  $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results/$paired_trimmed_read1_basename.bwa.sai \
			  $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results/$paired_trimmed_read2_basename.bwa.sai \
			  $paired_trimmed_read1 \
			  $paired_trimmed_read2 \
			  > $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results/$sample_folder\_bwa_aln_pe.sam
	
	cd $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results
	samtools view -bS $sample_folder\_bwa_aln_pe.sam > $sample_folder\_bwa_aln_pe.bam
	samtools sort -@ 8 -m 4G -o $sample_folder\_bwa_aln_pe.sorted.bam \
                     $sample_folder\_bwa_aln_pe.bam
	# index the big merged and sorted bam file
	samtools index $sample_folder\_bwa_aln_pe.sorted.bam 
	fi
	
done
	
