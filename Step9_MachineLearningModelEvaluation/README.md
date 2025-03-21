# Penetrance Prediction Model Validation

The machine learning penetrance prediction models used in the study are included in the Models directory. They are also compressed for easy upload to remote servers. There are three scripts for that produce the penetrance predictions in the validation dataset (in the case of the manuscript, All of Us). Note, the datasets for validation can be created using the CreateTraining_ValidationDatasets.py from Step8_MachineLearningModelFitting.

Step 1: Predict the penetrance scores for the pLoFs in the validation datasets (PredictFrameshiftPenetrance.py, PredictSpliceChangePenetrance.py, PredictStopGainPenetrance.py).

Step 2: Evaluate the machine learning model predictive performance on the disease expression measurements (MachineLearningModelEvaluation.py).


