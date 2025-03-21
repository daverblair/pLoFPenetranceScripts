import pandas as pd 
from tqdm import tqdm
import numpy as np
import pickle
from scipy import sparse


hpo_dx_table=pd.read_pickle('/Path/to/Clinical/Data/HPODx.pth')
all_hpos=sorted(list(set().union(*hpo_dx_table['HPO'].values)))
hpo_conversion=pd.Series(np.arange(len(all_hpos)),index=all_hpos)

data_vals=[]
row_inds=[]
col_inds=[]

subject_idx_counter=0
subject_ids=[]
for idx,row in hpo_dx_table.iterrows():
    hpo_dx_idx=[hpo_conversion.loc[x] for x in row['HPO']]]
    data_vals+=[1]*len(hpo_dx_idx)
    row_inds+=[subject_idx_counter]*len(hpo_dx_idx)
    col_inds+=hpo_dx_idx
    subject_idx_counter+=1
    subject_ids+=[str(idx)]
sparse_data_array=sparse.csr_matrix((data_vals,(row_inds,col_inds)),shape=(len(subject_ids),len(hpo_conversion)))
subject_ids=pd.Series(np.arange(len(hpo_dx_table)),index=subject_ids)

data_storage_dict={'SubjectIndex':subject_ids,'HPOColumns':hpo_conversion,'SparseSymptomMatrix':sparse_data_array}
with open('/Path/to/Clinical/Data/HPODx_SparseCSR.pth','wb') as f:
    pickle.dump(data_storage_dict,f)


