# Build Carrier and Non-Carrier Datasets

This script in this directory (BuildCarrier_NonCarrier_Datasets.py) builds a unique table for each disease, which has one entry per pLoF-variant carrier pair (carriers with multiple variants are dropped from the analysis). Each entry includes all of the variant-specific annotations required penetrance estimation/prediction. In addition, it constructs a list-like data structure containing all the non-carriers and other rare variant carriers. All of this information is stored in a pickled python object for downstream use. This script requires several other data files, described as follows:

1) Sample Count Files: These files contain the number of rare variants carried by each biobank subject, generated using plink2. See 'compute_sample_counts.sh' for an example script for computing these counts in the UKBB RAP. These counts are used to identify the non-carriers. Note, a different strategy was used for AoU based on hail.
2) The table of variant carriers created in 'Step1_LoF_Identification/IdentifyLoFCarriers.py'

3) The pLoF variant genomic features/scores created by the scripts in Step2_LoF_Annotation directory

4) A list of biobank subjects that pass exome level QC (see main methods)

5) A tab-delimited text file/table that contains the clinical data coverage statistics for the biobank. This requires querying the OMOP tables specific to each biobank.



