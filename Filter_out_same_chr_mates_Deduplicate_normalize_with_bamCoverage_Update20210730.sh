#!/bin/bash
# conda activate NGS_Py2.7 # installed in this env
 

top_level_folder="/data/guang/Experiment_20210312/raw_data/Zht8C_CT"

cd $top_level_folder
# sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"

sample_folder_names="$(ls -l | grep "^d" | grep "Z8" | awk '{print $NF}')" # make sure only folders with sample data will be selected

for sample_folder in $sample_folder_names
do
	
	#get into sub_directory of sample_folder which contains the bam file 
	cd $top_level_folder/$sample_folder/$sample_folder\_bwa_alignment_results
	
	# filter out reads with both mates from same chromosome using the 7th column "=" label
	grep '^@SQ' $sample_folder\_bwa_aln_pe.sam > $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped0.sam
	grep '^@PG' $sample_folder\_bwa_aln_pe.sam >> $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped0.sam
	# uniq mapping; 2 pair-mates from same chromosome; mapq >5
	# grep "XT:A:U" $sample_folder\_bwa_aln_pe.sam | awk '($7=="=" && $5>5)'  >> $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sam
	# select properly-mapped pair-mates(0x2) and slelect unique mapping (XT:A:U) and make sure pair-mates align to same chromosome ($7=="=")
	samtools view -hf 0x2 $sample_folder\_bwa_aln_pe.sam | grep "XT:A:U" | awk '$7=="="' >> $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped0.sam
	# grep XT:A:U $sample_folder\_bwa_aln_pe.sam >> $sample_folder\_bwa_aln_pe.uniq.mapped.sam
	
	# filter out reads with two mates, again.(previous grep and awk may create reads without mates)
	grep '^@SQ' $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped0.sam > $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sam
	grep '^@PG' $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped0.sam >> $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sam
	samtools view -f 0x2 $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped0.sam >> $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sam
	
	# Transfer sam into bam and sort index
	samtools view -bS $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sam > $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.bam
	samtools sort -@ 8 -m 4G -o $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.bam \
                     $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.bam
	# index the big merged and sorted bam file
	samtools index $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.bam
	# remove duplicates
	java -jar -Xmx8g /home/guang/bio_softwares/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
	METRICS_FILE=$sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.bam.dup.txt \
	VALIDATION_STRINGENCY=LENIENT \
	INPUT=$sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.bam \
	OUTPUT=$sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.bam
	# MarkDuplicates -REMOVE_DUPLICATES true -METRICS_FILE  Z8CT1_bwa_aln_pe.sorted.bam.dup.txt -INPUT Z8CT1_bwa_aln_pe.sorted.bam -OUTPUT Z8CT1_bwa_aln_pe.sorted.bam.deduplicated.bam
	
	# index deduplicated bam
	samtools sort -@ 8 -m 4G -o $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted.bam \
                     $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.bam
	samtools index $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted.bam
	
	# Normalize (deduplicated bam) with RPKM 
	bamCoverage -b $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted.bam \
	--normalizeUsing RPKM \
	--numberOfProcessors 8 \
	-o $sample_folder\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted_RPKM_normalized_noExtendReads.bw

done

