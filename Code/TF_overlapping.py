import argparse
from collections import defaultdict
import subprocess
import numpy as np
import re



parser=argparse.ArgumentParser()
parser.add_argument('-i',help='the input file',required=True)
parser.add_argument('-o',help='the ouput file',required=True)
parser.add_argument('-p', help = 'the path to the tf peak file')
parser.add_argument('-c', help = 'the col of chrom, enhancer_start, enhancer_end, promoter_start, promoter_end', type = int, nargs = '+', default = [0,1,2,4,5])
args=parser.parse_args()



def GetSlurmList(path):
   if path[-1] != '/': path += '/'
   return [path+i for i in subprocess.check_output(['ls',path]).decode('utf-8').split('\n') if re.search('ENCFF',i)]



def store_file(tf_file):
   tf= defaultdict(list)
   with open(tf_file) as f:
      for line in f:
         chrom, start, end = line.strip().split('\t')[:3]
         tf[chrom].append([int(start),int(end)])
   return Sorted(tf)


def check_mapping(store,regulatory_element,tf,tf_name):  ##peak here refers to tf
   for chrom in regulatory_element:
      if chrom not in tf: pass 
      else:
         pointer=0
         peak=tf[chrom]
         total_peak=len(peak)
         for element in regulatory_element[chrom]:
            start, end =element
            if (end < peak[pointer][0]): pass
            elif (start>peak[pointer][1]):
               if (pointer==(total_peak-1)): pass
               else:
                  while (start>peak[pointer][1] and pointer< (total_peak-1)): pointer+=1
                  if (end<peak[pointer][0]): pass
                  elif (start >peak[pointer][1] ): pass
                  else: store[chrom][element].append(tf_name)  ###don't need to consider multiple tf peaks overlapping the same regulatory elements, as only record if this TF overlapping the RE or not
            else: store[chrom][element].append(tf_name)
   return store   

def Sorted(input_dict):
   output_dict = {}
   for a_key in input_dict: output_dict[a_key] = sorted(input_dict[a_key])
   return output_dict

def TFpairs(list1, list2):
   pairs = set()
   set1, set2 = set(list1), set(list2)
   for tf1 in set1:
      for tf2 in set2: pairs.add((tf1,tf2))
   return list(pairs)




##########main function

input_file, columns, output_file, path  = args.i, args.c, args.o, args.p

#input_file, columns, output_file, path  = '/Users/hbb/Dropbox/PPI_project/Results/interactions/GM12878/MikeSnyder_GM12878_Rad21_positive.txt', [0,1,2,4,5], 'a', '/Users/hbb/Dropbox/PPI_project/Results/TF_peak_files/K562/test'
enhancer, promoter = defaultdict(dict), defaultdict(dict)
with open(input_file) as f:
   for line in f:
      chrom, enhancer_start, enhancer_end, promoter_start, promoter_end = np.array(line.strip().split('\t'))[columns]
      enhancer[chrom][(int(enhancer_start),int(enhancer_end))] = []
      promoter[chrom][(int(promoter_start),int(promoter_end))] = []


enhancer_pos_list = Sorted(enhancer)
promoter_pos_list = Sorted(promoter)



tf_list = GetSlurmList(path)
for file_name in tf_list:
   tf=store_file(file_name)
   tf_name= file_name.split('/')[-1].split('_')[0]
   enhancer=check_mapping(enhancer,enhancer_pos_list,tf,tf_name)
   promoter=check_mapping(promoter,promoter_pos_list,tf,tf_name)


f, w = open(input_file), open(output_file,'w')
for line in f:
   chrom, enhancer_start, enhancer_end, promoter_start, promoter_end = np.array(line.strip().split('\t'))[columns]
   enhancer_promoter_pair = TFpairs(enhancer[chrom][(int(enhancer_start),int(enhancer_end))], promoter[chrom][(int(promoter_start),int(promoter_end))])
   if enhancer_promoter_pair == []:  pass
   else: 
      w.writelines('\t'.join([chrom, enhancer_start, enhancer_end, promoter_start, promoter_end]))
      for a_pair in enhancer_promoter_pair: w.writelines('\t'+a_pair[0]+','+a_pair[1])
      w.writelines('\n')

f.close()
w.close()



