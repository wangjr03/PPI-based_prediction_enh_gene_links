




import argparse
import pandas as pd
import numpy as np
import os
import scipy.cluster.hierarchy as spc
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from scipy import interp
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve,auc
from sklearn.model_selection import StratifiedKFold
import matplotlib.patches as patches
import numpy as np 
import pandas as pd
from scipy.stats import pearsonr
import re
import random
from sklearn.feature_selection import RFECV

def Parser():
   parser=argparse.ArgumentParser('')
   parser.add_argument('-c', help = 'the result of the hierarchical TF community detection')
   parser.add_argument('-t', help = 'the path to the folder containing TF ChIP-seq peaks')
   parser.add_argument('-fn', help = 'the file of the feature names')
   parser.add_argument('-fm', help = 'the file of the feature matrix')
   parser.add_argument('-p', help = 'the path to the postive samples')
   parser.add_argument('-n', help = 'the path to the negative samples')
   parser.add_argument('-l', help = 'the path to the label list of the training samples')
   parser.add_argument('-o', help = 'output path of the trained model')
   parser.add_argument('-s', help = 'suffix')
   return parser.parse_args()   

args = Parser()
print('successfully read parameters')
labels = pd.read_csv(args.c,sep="\t",header=None)
print(labels)
all_TF_name = list(labels[0])
labels_2 = list(labels[2])
labels = list(labels[1])
TF_name = all_TF_name
print('Got datasets')




#get clusters of TFs
cluster = {}

for i in range(len(TF_name)):
    
    cluster[TF_name[i]] = labels[i]
    

    
    
    
#for those not in the cluster, add them:
#binding_list = os.listdir("/mnt/gs18/scratch/users/wangha73/PPI_data/TF_peak_files_2/GM12878/")
binding_list = os.listdir(args.t)

binding_list = [f for f in binding_list if '.txt' in f]

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

binding_list = [f for f in binding_list if '.txt' in f]

TF_name = [f.split("_")[0] for f in binding_list]

tmp_group = np.max(labels_2)+1
for i in TF_name:
    if i not in cluster_2:
        cluster_2[i] = tmp_group
        tmp_group += 1

        
        
        
        

        





#modify features: cluster TFs based on clusters

feature_name = pd.read_csv(args.fn)
feature_name=feature_name.drop('Unnamed: 0',axis=1)

feature_matrix = pd.read_csv(args.fm)
feature_matrix=feature_matrix.drop('Unnamed: 0',axis=1)

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


        
        
        
        
        
        
        


mapping = pd.read_csv(args.c,sep="\t")

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




data_split = 30






pos_pair= pd.read_csv(args.p,sep="\t",usecols = [0,1,2,3,4],header=None)
neg_pair = pd.read_csv(args.n,sep="\t",usecols = [0,1,2,3,4],header=None)


pos_pair['labels'] = 1
neg_pair['labels'] = -1

all_pair = pd.concat([pos_pair,neg_pair],axis=0,ignore_index=True)
all_pair.index = list(range(0,all_pair.shape[0]))

all_pair = all_pair.sort_values(by=[0,1])

all_pair.index = [list(range(0,all_pair.shape[0]))]
print(all_pair)

all_enh_list = np.array(all_pair.iloc[:,[0,1,2]])
    
all_pro_list = np.array(all_pair.iloc[:,[0,3,4]])
    
enh_dic = {}
pro_dic = {}
    
for i in range(len(all_enh_list)):
        
    if tuple(all_enh_list[i]) not in enh_dic:
        enh_dic[tuple(all_enh_list[i])] = [tuple(all_pro_list[i])]
    else:
        enh_dic[tuple(all_enh_list[i])].append(tuple(all_pro_list[i]))
 
    if tuple(all_pro_list[i]) not in pro_dic:
        pro_dic[tuple(all_pro_list[i])] = [tuple(all_enh_list[i])]
    else:
        pro_dic[tuple(all_pro_list[i])].append(tuple(all_enh_list[i]))



#build a set with non-overlapping enhancers and genes
#first define a set of consequetive enhancers and then find all genes linking to these enhancers

u_enh_list = all_pair.iloc[:,[0,1,2]].drop_duplicates()
u_pro_list = all_pair.iloc[:,[0,3,4]].drop_duplicates()

u_enh_list.index = [list(range(0,u_enh_list.shape[0]))]

u_enh_list.columns = ['enhancer_chrom','enhancer_start','enhancer_end']

#enh_idx = [int(f) for f in (np.array(u_enh_list.iloc[:,1]/1e7))]
enh_d = []

for i in range(u_enh_list.shape[0]-1):
    if u_enh_list.iloc[i,0] != u_enh_list.iloc[i+1,0]:
        enh_d.append( float('inf') )
    else:
        enh_d.append( (u_enh_list.iloc[i+1,1] + u_enh_list.iloc[i+1,2])/2 -(u_enh_list.iloc[i,1] + u_enh_list.iloc[i,2])/2 )

        
block_id = [1]
tmp_id = 1
for i in enh_d:
    
    if i < 1e6:
        block_id.append(tmp_id)
    else:
        tmp_id += 1
        block_id.append(tmp_id)

        
enh_idx = block_id


#devide them into 5 groups with similar number of enhancers

count_dic = {}

for i in enh_idx:
    if i in count_dic:
        count_dic[i] += 1
    else:
        count_dic[i] = 1

def get_test_set(count_dic,enh_idx,all_enh_list,frac):
    
    rd_key = random.sample(count_dic.keys(),len(count_dic.keys())-1)
    
    count = 0
    
    used_key = []
    
    for i in rd_key:
        count += count_dic[i]
        if count > frac*len(u_enh_list):
            return(used_key)
        else:
            used_key.append(i)
    





#this part will discover which feature to merge, according to out of bag accuracy and model complexity
#model complecity is estimated based on GDF:
#Computing AIC for black-box models using Generalised Degrees of Freedom: a comparison with cross-validation
def perturb_data(y,frac):
    #y is label list
    #frac is fraction to be perturbed
    y = pd.DataFrame(y)
    sel_ele = y.sample(frac=frac)
    new_label = y.copy()
    for i in sel_ele.index:
        tmp = 0-new_label.iloc[i,0]
        new_label.iloc[i,0] = tmp
    
    sel_loc = list(sel_ele.index)
    return sel_loc,np.array(new_label)
    
            
    
    
    

def GDF(clf,n,x,y,frac):
    #clf is the defined model
    #n is number of perturbation
    #x is feature matrix
    #y is label matrix
    #frac is fraction of 1 to be perturbed, same fraction will be used of 0
    lm = clf.fit(x,y).predict_proba(x)
    gdf = []
    for i in range(n):
        pt_loc,pt_y = perturb_data(pd.DataFrame(y),frac)
        y_flatten = pt_y.flatten()
        pt_score = clf.fit(x,y_flatten).predict_proba(x)
        gdf_up = [pt_score[i][1] - lm[i][1] for i in pt_loc]
        gdf_down = [pt_y[i][0] - y[i][0] for i in pt_loc]
        gdf.append(np.sum(gdf_up)/(1+np.sum(gdf_down)))
    return np.mean(gdf)

def AIC(clf,x,y,frac,n):
    gdf = GDF(clf,n,x,y,frac)
    clf.fit(x,y).predict_proba(x)
    lm_score = clf.oob_score_
    return -2*lm_score+2*gdf+gdf*(gdf+1)/(x.shape[1]-gdf-1)



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


        

        
        
        
feature_list = new_feature_mat
label_list = pd.read_csv(args.l)
label_list = label_list.drop('Unnamed: 0',axis=1)
        
x, y  = np.array(new_feature_mat), np.array(label_list)           
#for each directed data, merge it and compare AIC
clf = RandomForestClassifier(n_estimators=80,oob_score=True)
init_aic = AIC(clf,x,y,0.01,10)
tmp_y = [i[0] for i in y]
for i in rf_dic:
    if len(rf_dic[i]) == 1:
        continue
    else:
        tmp_feature = new_feature_mat.copy()
        tmp_col = 1*(new_feature_mat[rf_dic[i]].sum(axis=1)>0)
        tmp_feature = tmp_feature.drop(rf_dic[i],axis=1)
        tmp_col_name = list(tmp_feature.columns)
        tmp_feature = pd.concat([tmp_feature,tmp_col],axis=1)
        tmp_col_name.append(i)
        tmp_feature.columns = tmp_col_name
        aic = AIC(clf,np.array(tmp_feature),y,0.01,10)
        print(aic < init_aic)
        if aic < init_aic:
            new_feature_mat = tmp_feature             
            #init_aic = aic
        
        

test_list = []
train_list = []


for i in range(30):
    test_enh = get_test_set(count_dic,enh_idx,all_enh_list,0.2)
    test_enh = [f for f in range(len(enh_idx)) if enh_idx[f] in test_enh]
    test_enh = np.array(u_enh_list.iloc[test_enh,:])
    test_enh = [tuple(f) for f in test_enh]
    train, test = [],[]
    all_pair = pd.concat([pos_pair,neg_pair])
    all_pair.index = list(range(all_pair.shape[0]))
    #all_pair = all_pair_bp
    #all_pair = all_pair[rm_idx]
    for j in range(all_pair.shape[0]):
        if tuple(all_pair.iloc[j,0:3]) in test_enh:
            test.append(j)
        else:
            train.append(j)
    test_list.append(test)
    train_list.append(train)

aucs = []
tprs, aucs = [], []
importance = []
mean_fpr = np.linspace(0,1,100)
fi = []
feature_list = new_feature_mat
all_tpr = []
all_fpr = []
all_score = []
all_test_id=[]
for i in range(len(test_list)):  
    clf = RandomForestClassifier(n_estimators=80)
    mean_fpr = np.linspace(0,1,100)
    x, y  = np.array(feature_list), np.array(label_list)        
    prediction = clf.fit(x[train_list[i]],y[train_list[i]]).predict_proba(x[test_list[i]])    
    fi.append(clf.feature_importances_)
    fpr, tpr, t = roc_curve(y[test_list[i]], prediction[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    all_score.extend(prediction[:,1])
    all_test_id.extend(test_list[i])
    roc_auc = auc(fpr, tpr) 
    aucs.append(roc_auc)
    importance.append(clf.feature_importances_)

    
plt.plot([0,1],[0,1],linestyle = '--',lw = 2,color = 'black')
mean_tpr = np.mean(tprs, axis=0)
mean_auc = auc(mean_fpr, mean_tpr)
#OutputImportance(feature_space, importance, output_file)
plt.plot(mean_fpr, mean_tpr, color='red', label=r'Mean ROC (AUC = %0.2f )' % (mean_auc),lw=2, alpha=1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC')
plt.legend(loc='lower right',prop={'size': 6})
fi = pd.DataFrame([feature_list.columns,np.mean(importance,axis=0)]).transpose()

# train the model for later use.
import pickle as pk

final_model = RandomForestClassifier(n_estimators=80)
x, y  = np.array(feature_list), np.array(label_list)
final_model.fit(x,y)
pk.dump(final_model, open(args.o+"/Protect_model_"+args.s+".sav", 'wb'))
fi.to_csv(args.o+'/feature_importance_'+args.s+'.txt')
df = pd.DataFrame(list(zip(all_test_id, all_score)),columns =['row_id', 'score'])
df.to_csv(args.o+'/CV_scores_'+args.s+'.csv')
np.save(args.o+'TPR_'+args.s,mean_tpr)
np.save(args.o+'FPR_'+args.s,mean_fpr)
