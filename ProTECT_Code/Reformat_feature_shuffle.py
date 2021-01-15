
# coding: utf-8

# In[ ]:


import numpy as np
import os
import scipy.cluster.hierarchy as spc
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from ast import literal_eval        
import pandas as pd
import pickle as pk
import argparse



def Parser():
   parser=argparse.ArgumentParser('')
   parser.add_argument('-c', help = 'the result of the hierarchical TF community detection')
   parser.add_argument('-t', help = 'the path to the folder containing TF ChIP-seq peaks')
   parser.add_argument('-fn', help = 'the file of the feature names')
   parser.add_argument('m',help=  'path to the pre-trained model')
   parser.add_argument('-tn', help = 'the path to the file containing the feature name of the training data after feature engineering')
   parser.add_argument('-ifn', help = 'the path to the input feature name')
   parser.add_argument('-ifm', help = 'the path to the input feature matrix')
   parser.add_argument('-od', help= 'output path of the reformated features')
   parser.add_argument('-op', help='output path of the predicted propobilities to be postive')
   return parser.parse_args()   

args = Parser()



labels = pd.read_csv(args.c,sep="\t",header=None)
all_TF_name = list(labels[0])
labels_2 = list(labels[2])
labels = list(labels[1])
TF_name = all_TF_name
#get clusters of TFs
cluster = {}
for i in range(len(TF_name)):
    
    cluster[TF_name[i]] = labels[i]
    

    
    
    
#for those not in the cluster, add them:
binding_list = os.listdir(args.t)

binding_list = [f for f in binding_list if '.bed' in f]

TF_name = [f.split("_")[0] for f in binding_list]

tmp_group = np.max(labels)+1
for i in TF_name:
    if i not in cluster:
        cluster[i] = tmp_group
        tmp_group += 1

        
        
        
        
  


#repeat it for cluster 2
#get clusters of TFs
TF_name = all_TF_name
cluster_2 = {}

for i in range(len(TF_name)):
    
    cluster_2[TF_name[i]] = labels_2[i]
    

    
    
    
    
    
#for those not in the cluster, add them:
#binding_list = os.listdir("/mnt/gs18/scratch/users/wangha73/PPI_data/TF_peak_files_2/GM12878/")
binding_list = os.listdir(args.t)

binding_list = [f for f in binding_list if '.bed' in f]

TF_name = [f.split("_")[0] for f in binding_list]

tmp_group = np.max(labels_2)+1
for i in TF_name:
    if i not in cluster_2:
        cluster_2[i] = tmp_group
        tmp_group += 1


        

        
#modify features based on trained model        
sel_feature_space = pd.read_csv(args.fn)
sel_feature_vec = [i for i in sel_feature_space.iloc[:,0]]       

# In[ ]:


def process_chunk(feature_matrix,feature_name,train_feature_name,mapping):
    feature_matrix=feature_matrix.drop('Unnamed: 0',axis=1)
    feature_matrix.columns = list(feature_name['0'])
    for i in list(train_feature_name['0']):
        if i not in list(feature_name['0']):
            feature_matrix[i] = 0
    feature_matrix = feature_matrix[list(train_feature_name['0'])]
    feature_name = pd.DataFrame((feature_matrix.columns))
    TF_feature = [f for f in feature_name.iloc[:,0] if ',' in f]
    non_TF_feature =[f for f in feature_name.iloc[:,0] if ',' not in f]
    cluster_idx = []
    for i in TF_feature:
        TF1,TF2 = i.split(',')
        if TF1 in cluster and TF2 in cluster:
            cluster_idx.append([cluster[TF1],cluster[TF2]])
    #loop through it to generate a dictionary of columns
    col_dic = {}
    for i in range(len(cluster_idx)):
        if tuple(cluster_idx[i]) not in col_dic:
            col_dic[tuple(cluster_idx[i])] = [i]
        else:
            col_dic[tuple(cluster_idx[i])].append(i)
    new_feature_mat = pd.DataFrame()
    for i in col_dic:
        tmp_col = col_dic[i]
        new_feature = feature_matrix.iloc[:,tmp_col].sum(axis=1)
        new_feature_mat = pd.concat([new_feature_mat,1*(new_feature>0)],axis=1)
    #new_feature_mat = pd.concat([new_feature_mat,new_feature/new_feature.sum()],axis=1)
    new_feature_mat = pd.concat([new_feature_mat,feature_matrix.iloc[:,len(TF_feature):]],axis=1)
    a = list(col_dic.keys())
    a.extend(non_TF_feature)
    new_feature_mat.columns = a
    sel_col = []
    for i in a:
        if i in ['gene_activity', 'enhancer_activity', 'enhancer-gene_core', 'distance']:
            sel_col.append(i)
        elif i[0] == i[1]:
            sel_col.append(i)
    m_dic = {}
    for i in range(mapping.shape[0]):
        m_dic[mapping.iloc[i,0]] = mapping.iloc[i,1]
    within_cluster_feature = pd.DataFrame()
    for i in new_feature_mat.columns:
        if i[0] in m_dic and i[1] in m_dic and m_dic[i[0]] == m_dic[i[1]]:
             within_cluster_feature = pd.concat([within_cluster_feature,new_feature_mat[i]],axis=1)
    #within_cluster_feature = pd.concat([within_cluster_feature,new_feature_mat[['gene_activity','enhancer_activity','enhancer-gene_core','distance']]],axis=1)            
    cluster_idx = []
    for i in TF_feature:
        TF1,TF2 = i.split(',')
        if TF1 in cluster_2 and TF2 in cluster_2:
            cluster_idx.append([cluster_2[TF1],cluster_2[TF2]])
    #loop through it to generate a dictionary of columns
    col_dic = {}
    for i in range(len(cluster_idx)):
        if tuple(cluster_idx[i]) not in col_dic:
            col_dic[tuple(cluster_idx[i])] = [i]
        else:
            col_dic[tuple(cluster_idx[i])].append(i)
    new_feature_mat = pd.DataFrame()
    for i in col_dic:
        tmp_col = col_dic[i]
        new_feature = feature_matrix.iloc[:,tmp_col].sum(axis=1)
        new_feature_mat = pd.concat([new_feature_mat,1*(new_feature>0)],axis=1)
    #new_feature_mat = pd.concat([new_feature_mat,new_feature/new_feature.sum()],axis=1)
    new_feature_mat = pd.concat([new_feature_mat,feature_matrix.iloc[:,len(TF_feature):]],axis=1)
    a = list(col_dic.keys())
    a = [i+(1,) for i in a]
    a.extend(non_TF_feature)
    new_feature_mat.columns = a
    sel_col = []
    for i in a:
        if i[0] == i[1]:
            sel_col.append(i)
    new_feature_mat = new_feature_mat.drop(sel_col,axis=1)
    new_feature_mat = pd.concat([within_cluster_feature,new_feature_mat],axis=1)
    #merge reversed feature
    rf_dic = {}
    pool = list(new_feature_mat.columns)
    k=0
    assigned = []
    for i in new_feature_mat.columns[0:len(new_feature_mat.columns)-4]:
        tmp = (i[1],i[0])
        if i not in assigned:
            rf_dic[k] = [i]
            assigned.append(i)
            for j in new_feature_mat.columns:
                if (j[0],j[1]) == tmp and j not in assigned:
                    rf_dic[k].append(j)
                    assigned.append(j)
            k += 1
    #generate TF cluster features
    reformed_feature = []
    current_feature = [str(i) for i in list(new_feature_mat.columns)]
    assigned_feature = []
    for i in sel_feature_vec:
        if i in current_feature:
            loc = [j for j in range(len(current_feature)) if current_feature[j] == str(i)]
            assigned_feature.append(i)
            reformed_feature.append(list(new_feature_mat.iloc[:,loc[0]]))
    # generating merged features
    unassigned_feature = [i for i in sel_feature_vec if i not in assigned_feature]
    for i in unassigned_feature:
        tmp_col = rf_dic[int(i)]
        new_feature = 1*new_feature_mat[tmp_col].sum(axis=1)
        reformed_feature.append(list(new_feature))
    reformed_feature = pd.DataFrame(reformed_feature).transpose()
#convert columns name
    l = []
    for i in sel_feature_vec:
        if i not in ['gene_activity', 'enhancer_activity', 'enhancer-gene_core', 'distance']:
            l.append(literal_eval(i))
        else:
            l.append(i)
    reformed_feature.columns = l
    
    #shuffle features
    for i in reformed_feature.columns:
        reformed_feature[i] = reformed_feature[i].sample(frac=1).values
    
    return reformed_feature

    


loaded_model = pk.load(open(args.m, 'rb'))
train_feature_name = pd.read_csv(args.tn)
train_feature_name=train_feature_name.drop('Unnamed: 0',axis=1)
all_prob = []
mapping = pd.read_csv(args.c,sep="\t")
imp_feature_th = np.sort(loaded_model.feature_importances_)[-int(0.05*len(loaded_model.feature_importances_))]
l  = [i for i in range(len(loaded_model.feature_importances_)) if loaded_model.feature_importances_[i] > imp_feature_th]
all_flag = []
k = 1
def worker(suffix):
    feature_name = pd.read_csv(args.ifn)
    feature_name=feature_name.drop('Unnamed: 0',axis=1)
    feature_matrix_iter = pd.read_csv(args.ifm,chunksize=10000)
    k=1
    for chunk in feature_matrix_iter:
        processed_data = process_chunk(chunk,feature_name,train_feature_name,mapping)
        tmp_prob = loaded_model.predict_proba(processed_data)[:,1]
        #all_prob.extend(tmp_prob)
        tmp_prob = pd.DataFrame(tmp_prob)
        tmp_prob.to_csv(args.op,mode='a',index=None,header=None)
        #imp_f_name = processed_data.columns[l]
        #imp_f_name = imp_f_name.drop(['gene_activity','enhancer_activity', 'enhancer-gene_core','distance'])
        #flag = processed_data[imp_f_name].sum(axis=1)>0
        #all_flag.extend(list(flag))
        print(k)
        processed_data.to_csv(args.od,mode='a',index=None)
        k += 1

        
        
        
        
        
        
        
idx = list(range(0,28))
#idx = [0]
idy = []
for i in idx:
    if i < 10:
        idy.append('0'+str(i))
    else:
        idy.append(str(i))
        
for i in idy:
    worker(i)
    

