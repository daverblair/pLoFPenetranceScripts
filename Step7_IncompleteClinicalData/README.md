# Correcting for Low Data Coverage

These scripts repeat the symptom-driven penetance analysis, but first remove the subjects with low clinical data coverage.

Step 1-IdentifyLowCoverageSubjects.py: Uses logistic regression to develop a filter that removes subjects with very low clinical data coverage.


Step 2-pLoF_Penetrance_CombinedSymptomsDx_EHRFlaggedRemoved.py: Disease-specific pLoF-penetrance estimates are estimated from the global disease expression dataset, after removing subjects with low clinical data coverage.


