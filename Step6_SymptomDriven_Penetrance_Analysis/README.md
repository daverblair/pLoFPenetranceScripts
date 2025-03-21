# Symptom-Driven Penetrance Analysis

These scripts generate the results for our symptom-driven analyses. They are briefly described below:

Step 1-Compute_PheRS.py: Builds a matrix of Phenotype Risk Scores (PheRS). Each column corresponds to a disease, and each row represents a unique biobank subject.

Step 2-BrunnerMunzel_Analysis.py: This script compares the distribution of PheRS's observed in the carriers to those scores seen in the non-carriers.

Step 3-Estimtate_SympomDriven_Penetrance_Models.py: This script estimates the symptom-driven expression model and produces carrier-level expression probability estimates for each disease.

Step 4-Build_Global_Disease_Expression_Table.py: This script takes the results from Step 3 and builds a global table of disease expression measurments. It includes diagnoses, symptom-driven estimates, and their maximum. It also finds the threshold for symptom-driven expression score binarization (based on diagnoses) and binarizes the symptoms accordingly. 

Step 5-pLoF_Penetrance_Symptoms_Dx_Combined.py: Disease-specific pLoF-penetrance estimates are estimated from the global disease expression dataset. 


