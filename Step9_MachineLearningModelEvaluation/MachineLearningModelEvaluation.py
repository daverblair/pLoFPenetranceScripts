#!/usr/bin/env python
# coding: utf-8
import pandas as pd 
import numpy as np
from sklearn.utils import resample
import pickle

#This script analyzes the performance of the machine learning model using the measurements described in the main text. It requires 2 sets of input files: 1) the datasets of pLoF features used for penetrance prediction, and 2) the scores produced by the prediction models (see PredictFrameshiftPenetrance.py, PredictSpliceChangePenetrance.py, PredictStopGainPenetrance.py). It produces are the precision-recall curves and performance statistics displayed in Figure 5. 

stop_gain_annot_table = pd.read_pickle('/Path/to/ML/Model/Data/StopGainExpressionData_ValidationSet.pth')
frameshift_annot_table = pd.read_pickle('/Path/to/ML/Model/Data/FrameshiftExpressionData_ValidationSet.pth')
splice_annot_table = pd.read_pickle('/Path/to/ML/Model/Data/SpliceChangeExpressionData_ValidationSet.pth')

stop_gain_penetrance_scores = pd.read_pickle('/Path/to/ML/Model/Results/StopGainPenetranceScores_ValidationSet.pth')
frameshift_penetrance_scores = pd.read_pickle('/Path/to/ML/Model/Results/FrameshiftPenetranceScores_ValidationSet.pth')
splice_penetrance_scores = pd.read_pickle('/Path/to/ML/Model/Results/SpliceChangePenetranceScores_ValidationSet.pth')


def _custom_prec_recall_curve(penetrance_probs,variant_scores):
    desc_score_indices = np.argsort(variant_scores.values, kind="mergesort")[::-1]
    y_score=variant_scores.values[desc_score_indices]
    y_true=penetrance_probs.values[desc_score_indices]

    distinct_value_indices = np.where(np.diff(y_score))[0]
    threshold_idxs = np.r_[distinct_value_indices, y_true.size - 1]
    tps = np.cumsum(y_true)[threshold_idxs]
    fps = np.cumsum(1.0-y_true)[threshold_idxs]

    ps = tps + fps
    precision = np.zeros_like(tps)
    np.divide(tps, ps, out=precision, where=(ps != 0))
    recall = tps / tps[-1]
    sl = slice(None, None, -1)
    return np.hstack((precision[sl], 1)), np.hstack((recall[sl], 0)), y_score[threshold_idxs][sl]

def _custom_avg_precision_score(penetrance_probs,variant_scores):
    precision, recall, _ = _custom_prec_recall_curve(penetrance_probs,variant_scores)
    return max(0.0, -np.sum(np.diff(recall) * np.array(precision)[:-1]))


def _ppv_s_stats(observations,predictions):
    np.seterr(divide='ignore', invalid='ignore')
    p,r,t = _custom_prec_recall_curve(observations,predictions)
    f1 = 2*p*r/(p+r)
    max_f1 = np.nanmax(f1)
    loc_max_f1 = np.nanargmax(f1)
    thresh_max_f1 = t[loc_max_f1]
    ppv_max_f1 = p[loc_max_f1]
    sensitivity_max_f1 = r[loc_max_f1]
    avg_p_score = _custom_avg_precision_score(observations,predictions)
    return {'MaxF1':max_f1,'Score_MaxF1':thresh_max_f1,'PPV_MaxF1':ppv_max_f1,'Sens_MaxF1':sensitivity_max_f1,'Avg_PPV_Score':avg_p_score}
    
def _bootstrap_ppv_s_stats(observations,predictions,n_resamples=1000,ci=0.95):
    np.seterr(divide='ignore', invalid='ignore')
    running_means={'MaxF1':[],'Score_MaxF1':[],'PPV_MaxF1':[],'Sens_MaxF1':[],'Avg_PPV_Score':[]}
    for i in range(n_resamples):
        new_obs,new_pred = resample(observations,predictions)
        resampled_stats = _ppv_s_stats(new_obs,new_pred)
        for stat in running_means.keys():
            running_means[stat]+=[resampled_stats[stat]]
    for stat in running_means.keys():
        running_means[stat]=np.sort(np.array(running_means[stat]))
        
    low = (1.0-ci)/2
    high = 1.0-(1.0-ci)/2
    
    output_cis=dict([(x,None) for x in running_means.keys()])
    for key in output_cis.keys():
        output_cis[key] = (running_means[key][int(np.ceil(n_resamples*low))], running_means[key][int(np.ceil(n_resamples*high))])
    return output_cis

def _stat_perf_better_than_random(observations,predictions,n_resamples=10000,alternative='greater'):
    output_table={'BaselineAvgPrec':pd.NA,'ScoreAvgPrec':pd.NA,'PValue':pd.NA}
    baseline = observations.mean()
    score_avg_prec= _custom_avg_precision_score(observations,predictions)
    output_table['BaselineAvgPrec']=baseline
    output_table['ScoreAvgPrec']=score_avg_prec
    
    obs_diff = score_avg_prec-baseline
    num_at_least_as_extreme=0
    for i in range(n_resamples):
        new_obs = resample(observations,replace=False)
        new_baseline = new_obs.mean()
        new_score_avg_prec =  _custom_avg_precision_score(new_obs,predictions)
        if alternative=='greater':
            if (new_score_avg_prec-new_baseline)>=obs_diff:
                num_at_least_as_extreme+=1
        elif alternative=='less':
            if (new_score_avg_prec-new_baseline)<=obs_diff:
                num_at_least_as_extreme+=1
        else:
            if np.abs(new_score_avg_prec-new_baseline)>=np.abs(obs_diff):
                num_at_least_as_extreme+=1
    pval=num_at_least_as_extreme/n_resamples
    output_table['PValue']=pval
    return output_table

def _max_f1(observations,predictions):
    np.seterr(divide='ignore', invalid='ignore')
    p,r,t = _custom_prec_recall_curve(observations,predictions)
    f1 = 2*p*r/(p+r)
    max_f1 = np.nanmax(f1)
    return max_f1

def _recall_max_f1(observations,predictions):
    np.seterr(divide='ignore', invalid='ignore')
    p,r,t = _custom_prec_recall_curve(observations,predictions)
    f1 = 2*p*r/(p+r)
    max_f1_loc = np.nanargmax(f1)
    return r[max_f1_loc]

def _precision_max_f1(observations,predictions):
    np.seterr(divide='ignore', invalid='ignore')
    p,r,t = _custom_prec_recall_curve(observations,predictions)
    f1 = 2*p*r/(p+r)
    max_f1_loc = np.nanargmax(f1)
    return p[max_f1_loc]

def _bootstrap_compare_predictions(observations,predictions_1,predictions_2,stat_func,num_resamples=10000,alternative='greater'):
    output_table={'ObsDiff':pd.NA,'PValue':pd.NA}
    output_table['ObsDiff']=stat_func(observations,predictions_1)-stat_func(observations,predictions_2)
    
    num_failures=0
    for i in range(num_resamples):
        resampled_obs,resampled_pred_1,resampled_pred_2 = resample(observations,predictions_1,predictions_2)
        resampled_stat = stat_func(resampled_obs,resampled_pred_1)-stat_func(resampled_obs,resampled_pred_2)
        if alternative=='greater':
            if resampled_stat<0.0:
                num_failures+=1
        elif alternative=='less':
            if resampled_stat>0.0:
                num_failures+=1
 
    output_table['PValue']=num_failures/num_resamples
    return output_table 


shared_columns = ['SUBJECT_ID', 'DIS','ONSET','LOF_TYPE', 'SYMBOL','CLINVAR_ANNOT', 'CLINVAR_RATING','LOFTEE_CLASS','TX_ANNOT','PROB_EXPRESSION']
all_variant_types = pd.concat([stop_gain_annot_table[shared_columns],frameshift_annot_table[shared_columns],splice_annot_table[shared_columns]],axis=0)



all_penetrance_scores = pd.concat([stop_gain_penetrance__scores,frameshift_penetrance__scores,splice_penetrance__scores],axis=0)
all_variant_types=pd.concat([all_variant_types,all_penetrance_scores],axis=1)
all_variant_types=all_variant_types.loc[pd.isna(all_variant_types.RF_Score)==False]




### baseline filters ####
mane_loftee_filter = all_variant_types.index[(all_variant_types.TX_ANNOT=='MANE_SELECT')*(all_variant_types.LOFTEE_CLASS=='HC')]
mane_loftee_scores=pd.Series(np.zeros(all_variant_types.shape[0]),index=all_variant_types.index)
mane_loftee_scores[mane_loftee_filter]=1.0

plp_filter = all_variant_types.index[(all_variant_types.CLINVAR_ANNOT=='Pathogenic/Likely_Pathogenic')]
plp_scores=pd.Series(np.zeros(all_variant_types.shape[0]),index=all_variant_types.index)
plp_scores[plp_filter]=1.0

#### disease annotations ####

childhood_diseases = all_variant_types.index[all_variant_types.ONSET =='Childhood']
young_adult_diseases = all_variant_types.index[all_variant_types.ONSET =='YoungAdulthood']
adult_diseases = all_variant_types.index[all_variant_types.ONSET =='Adulthood']
stop_gains = all_variant_types.index[all_variant_types.LOF_TYPE =='STOP_GAINED']
frameshifts = all_variant_types.index[all_variant_types.LOF_TYPE =='FRAMESHIFT']
splice_variants = all_variant_types.index[all_variant_types.LOF_TYPE =='SPLICE_CHANGE']



# performance table
perf_stats_table = pd.DataFrame([],columns = ['MaxF1','Score_MaxF1','PPV_MaxF1','Sens_MaxF1','Avg_PPV_Score','Randomization_PValue'],index=['MANE_LoFTee','P_LP','ML_All','ML_StopGain','ML_Frameshift','ML_Splice','ChildhoodOnset','YoungAdultOnset','AdultOnset','ML_PLP_Only'])



# first baseline stats
mane_loftee_stats = _ppv_s_stats(all_variant_types.PROB_EXPRESSION,mane_loftee_scores)
mane_loftee_stats_95_ci = _bootstrap_ppv_s_stats(all_variant_types.PROB_EXPRESSION,mane_loftee_scores)
mane_loftee_pval = _stat_perf_better_than_random(all_variant_types.PROB_EXPRESSION,mane_loftee_scores)
for key in mane_loftee_stats.keys():
    perf_stats_table.at['MANE_LoFTee',key]={'Mean':mane_loftee_stats[key],'95_CI':mane_loftee_stats_95_ci[key]}
perf_stats_table.at['MANE_LoFTee','Randomization_PValue']=mane_loftee_pval


plp_stats = _ppv_s_stats(all_variant_types.PROB_EXPRESSION,plp_scores)
plp_stats_95_ci = _bootstrap_ppv_s_stats(all_variant_types.PROB_EXPRESSION,plp_scores)
plp_pval = _stat_perf_better_than_random(all_variant_types.PROB_EXPRESSION,plp_scores)
for key in plp_stats.keys():
    perf_stats_table.at['P_LP',key]={'Mean':plp_stats[key],'95_CI':plp_stats_95_ci[key]}
perf_stats_table.at['P_LP','Randomization_PValue']=plp_pval



ml_all_stats = _ppv_s_stats(all_variant_types.PROB_EXPRESSION,all_variant_types.RF_Score)
ml_all_95_ci = _bootstrap_ppv_s_stats(all_variant_types.PROB_EXPRESSION,all_variant_types.RF_Score)
ml_all_pval = _stat_perf_better_than_random(all_variant_types.PROB_EXPRESSION,all_variant_types.RF_Score)
for key in ml_all_stats.keys():
    perf_stats_table.at['ML_All',key]={'Mean':ml_all_stats[key],'95_CI':ml_all_95_ci[key]}
perf_stats_table.at['ML_All','Randomization_PValue']=ml_all_pval

precision,recall,thresholds = _custom_prec_recall_curve(all_variant_types.PROB_EXPRESSION,all_variant_types.RF_Score)
output_dictionary={}
output_dictionary['PrecisionScores']=precision
output_dictionary['RecallScores']=recall
output_dictionary['Thresholds']=thresholds
with open('/Path/to/ML/Model/Results/ML_All_PrecRecallCurve.pth','wb') as f:
    pickle.dump(output_dictionary,f)



ml_plp_stats = _ppv_s_stats(all_variant_types.loc[plp_filter].PROB_EXPRESSION,all_variant_types.loc[plp_filter].RF_Score)
ml_plp_95_ci = _bootstrap_ppv_s_stats(all_variant_types.loc[plp_filter].PROB_EXPRESSION,all_variant_types.loc[plp_filter].RF_Score)
ml_plp_pval = _stat_perf_better_than_random(all_variant_types.loc[plp_filter].PROB_EXPRESSION,all_variant_types.loc[plp_filter].RF_Score)
for key in ml_plp_stats.keys():
    perf_stats_table.at['ML_PLP_Only',key]={'Mean':ml_plp_stats[key],'95_CI':ml_plp_95_ci[key]}
perf_stats_table.at['ML_PLP_Only','Randomization_PValue']=ml_plp_pval

precision,recall,thresholds = _custom_prec_recall_curve(all_variant_types.loc[plp_filter].PROB_EXPRESSION,all_variant_types.loc[plp_filter].RF_Score)
output_dictionary={}
output_dictionary['PrecisionScores']=precision
output_dictionary['RecallScores']=recall
output_dictionary['Thresholds']=thresholds
with open('/Path/to/ML/Model/Results/ML_PLP_Only_PrecRecallCurve.pth','wb') as f:
    pickle.dump(output_dictionary,f)



indices=['ML_StopGain','ML_Frameshift','ML_Splice','ChildhoodOnset','YoungAdultOnset','AdultOnset']
variant_filters=[stop_gains,frameshifts,splice_variants,childhood_diseases,young_adult_diseases,adult_diseases]
                 
for i in range(len(indices)):
                 
    ml_stats = _ppv_s_stats(all_variant_types.loc[variant_filters[i]].PROB_EXPRESSION,all_variant_types.loc[variant_filters[i]].RF_Score)
    ml_95_ci = _bootstrap_ppv_s_stats(all_variant_types.loc[variant_filters[i]].PROB_EXPRESSION,all_variant_types.loc[variant_filters[i]].RF_Score)
    ml_pval = _stat_perf_better_than_random(all_variant_types.loc[variant_filters[i]].PROB_EXPRESSION,all_variant_types.loc[variant_filters[i]].RF_Score)
    for key in ml_stats.keys():
        perf_stats_table.at[indices[i],key]={'Mean':ml_stats[key],'95_CI':ml_95_ci[key]}
    perf_stats_table.at[indices[i],'Randomization_PValue']=ml_pval

    precision,recall,thresholds = _custom_prec_recall_curve(all_variant_types.loc[variant_filters[i]].PROB_EXPRESSION,all_variant_types.loc[variant_filters[i]].RF_Score)
    output_dictionary={}
    output_dictionary['PrecisionScores']=precision
    output_dictionary['RecallScores']=recall
    output_dictionary['Thresholds']=thresholds
    with open('/Path/to/ML/Model/Results/{0:s}_PrecRecallCurve.pth'.format(indices[i]),'wb') as f:
        pickle.dump(output_dictionary,f)



perf_stats_table.to_pickle('/Path/to/ML/Model/Results/OverallPrecisionRecallStats.pth')
perf_stats_table.to_csv('/Path/to/ML/Model/Results/OverallPrecisionRecallStats.txt',sep='\t')



# comparison analysis
ml_vs_labels=pd.DataFrame([],columns=['Avg_Precision','Max_F1','Recall_Max_F1','Precision_Max_F1'],index=['MANE_LoFTee','P_LP'])
func_list = [_custom_avg_precision_score,_max_f1,_recall_max_f1,_precision_max_f1]
func_labels = ml_vs_labels.columns

for i in range(len(func_list)):
    stats = _bootstrap_compare_predictions(all_variant_types.PROB_EXPRESSION,all_variant_types.RF_Score,mane_loftee_scores,func_list[i])
    ml_vs_labels.at['MANE_LoFTee',func_labels[i]]=stats
    
for i in range(len(func_list)):
    stats = _bootstrap_compare_predictions(all_variant_types.PROB_EXPRESSION,all_variant_types.RF_Score,plp_scores,func_list[i])
    ml_vs_labels.at['P_LP',func_labels[i]]=stats



ml_vs_labels.to_pickle('/Path/to/ML/Model/Results/Data/ML_vs_BasicLabelsStats.pth')
ml_vs_labels.to_csv('/Path/to/ML/Model/Results/Data/ML_vs_BasicLabelsStats.txt',sep='\t')






