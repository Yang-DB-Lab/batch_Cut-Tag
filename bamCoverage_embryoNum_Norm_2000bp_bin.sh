#!/bin/bash
top_level_folder="/data/guang/Experiment_20210312/raw_data/Zht8C_CT"

# folder to keep the bigwig and bedgraph files
mkdir -p /data/guang/Experiment_20210312/raw_data/Zht8C_CT/Zht_CutTag_DataAnalysis20210802/EmbryoNumber_normed_bw_2000bp_bin

cd $top_level_folder
# sample_folder_names="$(ls -l | grep "^d" | awk '{print $NF}')"

sample_folder_names=($(ls -l | grep "^d" | grep "Z8" | awk '{print $NF}'))
# make sure only folders with sample data will be selected
# add quote outside to transfer the results into an array !!!

scaleFactors=(6.0000000 6.0000000 2.0000000 1.9867550 1.7857143 1.4285714 1.0526316 0.8333333)
# array of scale factors(floating numbers)

for i in {0..7} # 0-based indexing
do
    #get into sub_directory of sample_folder which contains the bam file 
	cd $top_level_folder/${sample_folder_names[i]}/${sample_folder_names[i]}\_bwa_alignment_results


  bamCoverage -b ${sample_folder_names[i]}\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted.bam \
	--scaleFactor ${scaleFactors[i]} \
	--numberOfProcessors 8 \
	--binSize 2000 \
	--outFileFormat bedgraph \
	-o ${sample_folder_names[i]}\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted_EmbryoNumberNormed_noExtendReads.2000bp_bin.bedgraph
	
	cp ${sample_folder_names[i]}\_bwa_aln_pe.same_chromosome_mates.uniq.mapped.sorted.deduplicated.sorted_EmbryoNumberNormed_noExtendReads.2000bp_bin.bedgraph /data/guang/Experiment_20210312/raw_data/Zht8C_CT/Zht_CutTag_DataAnalysis20210802/EmbryoNumber_normed_bw_2000bp_bin
done
