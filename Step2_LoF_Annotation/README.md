# pLoF Annotation

These scripts annotate pLoF variants with the features described in the Supplementary Information. They require the output of the scripts in Step1_LoF_Identification.

ParseCADDScores.py:  Builds a table of CADD scores using the output from CADD's online scoring tool. 

ParseSpliceAIResults.py: Builds a table of the genomic features for the splice change variants. Uses functions defined in the Auxillary_Functions/SequenceAnalysisFunctions.py script.

PredictNMDEscape.py: Annotates NMD escape features for stop-gain variants. Builds table with results

PredictNumAAsImpacted_NonSplice.py: Annotates stop gain and splice change variants with the fraction of amino acids impacted. Splice change variants are annotated with this information by ParseSpliceAIResults.py


