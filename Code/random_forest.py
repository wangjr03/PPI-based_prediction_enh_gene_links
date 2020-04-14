'''
-i: the input files, nargs = '+'

-newin: the new input data

-newout: the new output data

-o: the output image name 

-l: the label of the input files, 1 for positive set and -1 for negative set, nargs = '+', type = int

-n: the number of trees in RF

-k: the k-fold cross-validation, default = 5

-s: output feature list in a file, with the feature importance in each CV

-activity:  input the cell line key word here. add enhancer activity, gene activity and gene-enhancer correlation

-d: store true, if add the distance feature or not
'''

import argparse
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


def Parser():
   parser=argparse.ArgumentParser('')
   parser.add_argument('-i', help = 'the input files', nargs = '+')
   parser.add_argument('-l', help = 'the labels of the corresponding input files. 1 for positive set and -1 for negative set', nargs = '+', type = int)
   parser.add_argument('-s', help = 'output feature list in a file, with the feature importance in each CV')
   parser.add_argument('-n', help = 'the number of trees in RF', type = int)
   parser.add_argument('-k', help = 'the k-fold cross-validation', type = int, default = 5)
   parser.add_argument('-o', help = 'the output image name')
   parser.add_argument('-activity', help = 'input the cell line key word here. add enhancer activity, gene activity and gene-enhancer correlation')
   parser.add_argument('-d', help = 'if add the distance feature or not', action = 'store_true')
   parser.add_argument('-gene_annotation', help = 'gene_annotation_file', default = '/mnt/research/compbio/wanglab/data/Roadmap/genes/Gene_annotation.bed')
   parser.add_argument('-gene_name', help = 'gene_name_file', default = '/mnt/research/compbio/wanglab/data/Roadmap/genes/RPKM_all_gene_name_select_2.bed')
   parser.add_argument('-cell_name', help = 'cell_name_file', default = '/mnt/research/compbio/wanglab/data/Roadmap/Cell_types_RNAseq.bed')
   parser.add_argument('-gene_activity', help = 'gene_activity_file', default = '/mnt/home/huangbi4/PPI_project/new/results/Roadmap/gene_activity_quantile_normalization.bed')
   parser.add_argument('-enhancer_activity', help = 'enhancer_activity file', default = '/mnt/home/huangbi4/PPI_project/new/results/Roadmap/enhancer_activity_quantile_normalization.bed')
   parser.add_argument('-enhancer_pos', help = 'enhancer_pos file', default = '/mnt/research/compbio/wanglab/data/Roadmap/enhancers/enhancer_coords_select.bed')
   parser.add_argument('-v4', help = 'the forth version, only using activity feature', action = 'store_true')
   parser.add_argument('-newin', help = 'the new input data')
   parser.add_argument('-newout', help = 'the new output data')
   return parser.parse_args()   



def GetCellLine(cell_list, key_word):
   for i in cell_list:
      if re.search(key_word.upper(), i.upper()): return i



def GeneToEnsembl(input_file):
   gene_ensembl = {}
   with open(input_file) as f:
      for line in f:
         line = line.split('\t')
         if line[4] == '1': gene_ensembl[('chr'+line[1], int(line[2]))] = line[0]
         else: gene_ensembl[('chr'+line[1], int(line[3]))] = line[0]
   return gene_ensembl



def CalculateCorrelation(list1, list2):
   return pearsonr(list1,list2)[0]   


def StoreFile(input_file, feature_space, feature_list, label, label_list):
   with open(input_file) as f:
      for line in f:
         label_list.append(label)
         features = line.strip().split('\t')[5:]
         feature_space.extend(features)
         feature_list.append(features)
   return feature_space, feature_list, label_list



def OutputImportance(feature_space,importance, output_file):
   with open(output_file, 'w') as w:
      for i in range(len(feature_space)):  ###for each feature
         w.writelines(feature_space[i])
         for j in range(len(importance)):  ###the importance in each CV
            w.writelines('\t'+str(importance[j][i]))
         w.writelines('\n')
   return 0




def ModifyFeatureList(feature_space_index, feature_list):
   store = []
   for a_list in feature_list:
      res = [0 for i in range(len(feature_space_index))]
      for a_feature in a_list:
         index =  feature_space_index[a_feature]
         res[index] = 1
      store.append(res)
   return store



def NewFeatureList(new_data, feature_space_index):
   store = []
   with open(new_data) as f:
      for line in f:
         tf_list = line.strip().split('\t')[5:]
         res = [0 for i in range(len(feature_space_index))]
         for a_feature in tf_list:
            try: 
               index = feature_space_index[a_feature]
               res[index] = 1
            except: pass
         store.append(res) 
   return store


def StoreDistance(store,input_file):
   with open(input_file) as f:
      for line in f:
         enhancer_start, enhancer_end, promoter_start, promoter_end = np.array(line.strip().split('\t'))[1:5].astype(int)
         store.append([abs(int((enhancer_start+enhancer_end)/2) - int((promoter_start+promoter_end)/2))])
   return store



def StoreActivity(activity_feature, input_file, cell_line, gene_activity, gene_activity_list, gene_name_index, gene_ensembl, enhancer_activity, enhancer_activity_list, enhancer_pos_index):
   with open(input_file) as f:
      for line in f:
         line = np.array(line.strip().split('\t'))
         chrom, [enhancer_start, enhancer_end, promoter_start, promoter_end] = line[0], line[1:5].astype(int)
         tss = int((promoter_start + promoter_end)/2)
         a_gene_activity = gene_activity[cell_line][gene_name_index[gene_ensembl[(chrom,tss)]]]
         a_enhancer_activity = enhancer_activity[cell_line][enhancer_pos_index[(chrom, enhancer_start, enhancer_end)]]
         pearson_cor = pearsonr(gene_activity_list[gene_name_index[gene_ensembl[(chrom,tss)]]], enhancer_activity_list[enhancer_pos_index[(chrom, enhancer_start, enhancer_end)]])[0]
         activity_feature.append([a_gene_activity, a_enhancer_activity, pearson_cor])
   return activity_feature



def OutputPrediction(new_data, new_output, new_prediction):
   f, w, index = open(new_data), open(new_output,'w'), 0
   prediction_round = len(new_prediction)
   for line in f:
      if sum([1 for i in range(prediction_round) if new_prediction[i][index] >0.5]) == prediction_round: w.writelines(line)
      index += 1
   f.close(), w.close()
   



def main():
   args = Parser()           
   input_files, labels, n_trees, k_fold, output_file, new_data, new_output = args.i, args.l, args.n, args.k, args.s, args.newin, args.newout
   feature_space, feature_list, label_list, activity_feature, distance, new_activity_feature, new_distance = [], [] ,[], [], [], [], []
######store files
   for i in range(len(input_files)): feature_space, feature_list, label_list = StoreFile(input_files[i], feature_space, feature_list, labels[i], label_list)
   print('finish store files')
######store features
   feature_space = list(set(feature_space))
   feature_space_index = dict(zip(feature_space, range(len(feature_space))))
   print('Number of features:'+str(len(feature_space)))
   feature_list = ModifyFeatureList(feature_space_index, feature_list)
   print('finish modify features')
#######store the features for new data
   new_feature_list = NewFeatureList(new_data, feature_space_index)
   if args.activity:
      gene_annotation_file, gene_name_file, cell_name_file, gene_activity_file, enhancer_activity_file, enhancer_pos_file = args.gene_annotation, args.gene_name, args.cell_name, args.gene_activity, args.enhancer_activity, args.enhancer_pos
      feature_space.extend(['gene_activity','enhancer_activity','enhancer-gene_core'])
      gene_name_list = pd.read_csv(gene_name_file, header = None, sep = '\t')[0].tolist()
      cell_list = pd.read_csv(cell_name_file, sep = '\t')['Universal_Human_Reference'].tolist()
      gene_name_index = dict(zip(gene_name_list,range(len(gene_name_list))))
      cell_line = GetCellLine(cell_list, args.activity)
      gene_ensembl = GeneToEnsembl(gene_annotation_file)
      gene_activity = pd.read_csv(gene_activity_file, header = None, sep = '\t', names = cell_list)
      gene_activity_list =  gene_activity.values.tolist()
      df = pd.read_csv(enhancer_pos_file, header = None, sep = '\t')
      enhancer_pos_list = list(zip(df[0],df[1],df[2]))
      enhancer_pos_index = dict(zip(enhancer_pos_list,range(len(enhancer_pos_list))))
      enhancer_activity = pd.read_csv(enhancer_activity_file, header = None, sep = '\t', names = cell_list)
      enhancer_activity_list = enhancer_activity.values.tolist()
      for i in range(len(input_files)): activity_feature = StoreActivity(activity_feature, input_files[i], cell_line, gene_activity, gene_activity_list, gene_name_index, gene_ensembl, enhancer_activity, enhancer_activity_list, enhancer_pos_index)
      new_activity_feature = StoreActivity(new_activity_feature, new_data, cell_line, gene_activity, gene_activity_list, gene_name_index, gene_ensembl, enhancer_activity, enhancer_activity_list, enhancer_pos_index)
      if args.v4: 
         feature_list = activity_feature
         feature_space = ['gene_activity','enhancer_activity','enhancer-gene_core']
         new_feature_list = new_activity_feature
      else: 
         feature_list = [feature_list[i]+activity_feature[i] for i in range(len(activity_feature))]     
         new_feature_list = [new_feature_list[i]+new_activity_feature[i] for i in range(len(new_activity_feature))] 
   if args.d:
      feature_space.append('distance')
      for i in range(len(input_files)): distance = StoreDistance(distance,input_files[i])
      new_distance = StoreDistance(new_distance, new_data)
      feature_list = [feature_list[i]+distance[i] for i in range(len(distance))]  
      new_feature_list = [new_feature_list[i]+new_distance[i] for i in range(len(new_distance))]    
######random forests
   clf = RandomForestClassifier(n_estimators=n_trees)
   cv = StratifiedKFold(n_splits=k_fold, shuffle=True)
   tprs, aucs = [], []
   #importance=clf.feature_importances_
   importance = []
   mean_fpr = np.linspace(0,1,100)
   i = 1
   x, y  = np.array(feature_list), np.array(label_list)
   new_prediction = []
   print('begin fitting model')
   for train,test in cv.split(x, y):
      print('round:'+str(i))
      prediction = clf.fit(x[train],y[train]).predict_proba(x[test])
      new_prediction.append([i[0] for i in clf.fit(x[train],y[train]).predict_proba(new_feature_list)])
      fpr, tpr, t = roc_curve(y[test], prediction[:, 1])
      tprs.append(interp(mean_fpr, fpr, tpr))
      roc_auc = auc(fpr, tpr)
      aucs.append(roc_auc)
      importance.append(clf.feature_importances_)
      plt.plot(fpr, tpr, lw=2, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
      i= i+1
   plt.plot([0,1],[0,1],linestyle = '--',lw = 2,color = 'black')
   mean_tpr = np.mean(tprs, axis=0)
   mean_auc = auc(mean_fpr, mean_tpr)
   OutputImportance(feature_space, importance, output_file)
   plt.plot(mean_fpr, mean_tpr, color='blue', label=r'Mean ROC (AUC = %0.2f )' % (mean_auc),lw=2, alpha=1)
   plt.xlabel('False Positive Rate')
   plt.ylabel('True Positive Rate')
   plt.title('ROC')
   plt.legend(loc="lower right")
   plt.savefig(args.o)
   OutputPrediction(new_data, new_output, new_prediction)
   return 0




############
main()


