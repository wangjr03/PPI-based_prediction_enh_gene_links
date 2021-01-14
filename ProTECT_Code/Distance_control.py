
# coding: utf-8

# In[ ]:



def Parser():
   parser=argparse.ArgumentParser('')
   parser.add_argument('-p', help = 'the path to the postive samples')
   parser.add_argument('-n', help = 'the path to the negative samples')
   parser.add_argument('-a', help = 'the output path to the corrected positive samples')
   parser.add_argument('-b', help = 'the out path to the corrected negative samples')
   return parser.parse_args()   

args = Parser()


import pandas as pd

pos = pd.read_csv(args.p,sep='\t',header=None)
neg = pd.read_csv(args.n,sep='\t',header=None)
pos_d = abs((pos.iloc[:,3]+pos.iloc[:,4])/2-(pos.iloc[:,1]+pos.iloc[:,2])/2)
neg_d = abs((neg.iloc[:,3]+neg.iloc[:,4])/2-(neg.iloc[:,1]+neg.iloc[:,2])/2)
pos['distance'] = pos_d
neg['distance'] = neg_d
pos = pos.sort_values(['distance'])
neg = neg.sort_values(['distance'])
pos_bin = []
neg_bin = []
for i in range(pos.shape[0]):
    pos_bin.append( int(abs((pos.iloc[i,3]+pos.iloc[i,4])/2-(pos.iloc[i,1]+pos.iloc[i,2])/2)/2e4) )
for i in range(neg.shape[0]):
    neg_bin.append( int(abs((neg.iloc[i,3]+neg.iloc[i,4])/2-(neg.iloc[i,1]+neg.iloc[i,2])/2)/2e4) )

    
    
    
    
pos_dic = {}
neg_dic = {}
for i in range(len(pos_bin)):
    if pos_bin[i] <= 50:
        if pos_bin[i] not in pos_dic:
            pos_dic[pos_bin[i]] = [i]
        else:
            pos_dic[pos_bin[i]].append(i)

            
            
for i in range(len(neg_bin)):
    if neg_bin[i] <= 50:
        if neg_bin[i] not in neg_dic:
            neg_dic[neg_bin[i]] = [i]
        else:
            neg_dic[neg_bin[i]].append(i)

pos_id = []
neg_id = []

for i in neg_dic:
    if i not in pos_dic:
        continue
    else:
        m = len(pos_dic[i])
        n = len(neg_dic[i])
        pos_id.extend(pos_dic[i][0:min(n,m)])
        neg_id.extend(neg_dic[i][max(0,n-int(1*min(n,m))):])

        
        
        
        
pos_set = pos.iloc[pos_id,]
neg_set = neg.iloc[neg_id,]

pos_d = abs((pos_set.iloc[:,3]+pos_set.iloc[:,4])/2-(pos_set.iloc[:,1]+pos_set.iloc[:,2])/2)
neg_d = abs((neg_set.iloc[:,3]+neg_set.iloc[:,4])/2-(neg_set.iloc[:,1]+neg_set.iloc[:,2])/2)

pos_set.to_csv(args.a,sep="\t",header=None,index=None)    
neg_set.to_csv(args.b,sep="\t",header=None,index=None)    

