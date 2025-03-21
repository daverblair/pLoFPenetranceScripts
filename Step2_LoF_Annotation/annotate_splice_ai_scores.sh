# conda activate tensorflow


#This script was used to generate the SpliceAI scores using custom transcript annots (ensembl 110)

fasta_fila=/Path/to/hg38/hg38.fa

input_file=/Path/to/LoF/VCF/AllHaploLOFVariants_NoSamples_Sorted.vcf.gz
output_file=/Path/to/SpliceAI/Output/AllHaploLOFVariants_spliceAI_Scores_500bp.vcf.gz
annot_file=/Path/to/CustomSpliceAIAnnots.txt

spliceai -I ${input_file} -O ${output_file} -R ${fasta_fila} -A ${annot_file} -D 500