# Create Clinical Datasets

These scripts create the clinical datasets used in the analysis of pLoF penetrance. The analysis  originates from a text file that contains all of the OMOP concept IDs assigned to every biobank participant, along with the dates that they were assigned. The precise code for generating this text file varies across biobanks, but since  both now use the OMOP-CDM format, the text files can be generated using an SQL query to the condition_occurrence table:     

SELECT condition_occurrence.subject_id, condition_occurrence.condition_concept_id, condition_occurrence.condition_start_date FROM data_set.condition_occurrence 

After this query, the results can be dumped into a text file, which can then be sorted by subject_id.


After creating this text file, all of the clinical datasets required for the analysis are created by running the following scripts:

Step 1-BuildOMOPCodeMatrix.py: This script takes the output from the text file and stores it as a compact, custom sparse data structure. 

Step 2-BuildDiseaseDxTables.py: From this sparse data structure, we create tables of haploinsufficient disease diagnoses for every biobank participant (one table for each disease).

Step 3-ConvertOMOP_to_HPO.py: The OMOP concept IDs are converted to HPO symptoms (if an aligned symptom exists).

Step 4-ConstructSparseHPOMatrices.py: The sparse data structure from the previous step is further converted into a scipy.sparse_csr matrix for additional analyses.


        


