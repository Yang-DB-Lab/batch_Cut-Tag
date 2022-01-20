#!/bin/bash
top_level_folder="/data/guang/Experiment_20210312/raw_data/Zht8C_CT"

cd $top_level_folder

for sample_folder in Z8CT{1..8}
do
	
	#get into sub_directory of sample_folder which contains the bam file 
	cd $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results
	
	samtools view $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted.bam \
	| cut -f9 | awk '$1>0' > $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted.bam.fragment.sizes.txt

done
