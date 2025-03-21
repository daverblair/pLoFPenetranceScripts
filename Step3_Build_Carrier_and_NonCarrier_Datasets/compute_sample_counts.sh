#!/usr/bin/env bash


#This script computes the variant counts for every subject in the UKBB using plink2 via the RAP. Note, you must first generate pgen variant files that are restricted to biallelic variants with an AF < 0.1% using standard tools. These are stored in a directory within the project called '/SequenceData/GenePGENS_RareOnlyBiallelic/' within the UKBB project directory


current_time=$(date +%s)
log_file_name=dnanexus_sample_counts_${current_time}.log
PROJECT=XXX

dx login --noprojects
dx select ${PROJECT} >> ${log_file_name} 2>&1


input_directory='/SequenceData/GenePGENS_RareOnlyBiallelic/'
output_directory=/SequenceData/SampleCounts/
gene_list="Path/to/Auxiullary_Data/GeneList.txt"

while IFS=$'\t', read -r SYMBOL; do
	echo "Performing sample counts for ${SYMBOL}..."

	FILE="${SYMBOL}_PlinkTmp"
	OUTPUT="${SYMBOL}_PerSampleCounts"
	command=" plink2 --pfile ${FILE} --sample-counts --out ${OUTPUT}"

	dx run swiss-army-knife -iin="${input_directory}${FILE}.pgen" -iin="${input_directory}${FILE}.psam" -iin="${input_directory}${FILE}.pvar" -imount_inputs=true -icmd="${command}" --tag="${SYMBOL}_PerSampleCounts" --destination=${output_directory} --instance-type="mem1_ssd1_v2_x8" >> ${log_file_name} 2>&1
	sleep 1

done < ${gene_list}

dx logout >> ${log_file_name} 2>&1


