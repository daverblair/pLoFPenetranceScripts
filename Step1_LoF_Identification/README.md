# Identify pLoFs and Carriers

These scripts assume that all biobank variants within the haploinsufficient disease genes of interest have been identified using standard tools. The precise tools required will vary according to the biobank (ex: bcftools in the UKBB vs hail/bcftools in All of Us) and researcher preferences. The variant data for this study was stored as follows. The sample-level information for each gene was removed the raw data, and all rare variants were distributed to gene-specific VCF files for annotation using VEP (ex: SYMBOL_Exome_Annotated_NoSamples.vcf.gz). A single, global VCF file with sample-level data was then constructed by isolating only pLoFs from the master VCF/hail table for the biobank (ex: AllHaploLOFVariants_wSamples_Sorted.vcf.gz).



Step 1-Identify pLoFs (IdentifyLoFs.py): This script identifies the pLoFs from a set of gene-specific VCF files that were stripped of sample level information (enabling them to be downloaded from biobank servers and analyzed locally). There are three outputs of this script: a table of all identified pLoFs along with basic annotation information (AllHaploLOFVariants_NoScores.txt), a list of all genes with pLoFs (GeneList_wLoF.txt), and gene-specific files containing the coordinates of all identified pLoFs (SYMBOL_HaploLOFVariants_TargetCoords.txt). The latter files were used in conjunction with bcftools/hail to isolate only the pLoFs of interest from the global variant data with sample information (i.e. build AllHaploLOFVariants_wSamples_Sorted.vcf.gz)

Step 2-Identify pLoF Carriers (IdentifyLoFCarriers.py): Using the global VCF file containing sample-level information for all pLoF variants identified in Step 1, this script pulls out all potential carriers, subjects the variants/calls to quality-control filtering, and stores all this information in a table for downstream analyes (AllHaploLOFVariants_FilteredCarriers.txt). Note, this script is for the UKBB specifically. Slight modifications would be required (see comments) for it to work in All of Us.


