'''
-i: the input file

-o: the output file

-s: the PPI score threshold, type = int, full_score = 1000

-p: the ppi file, default = '/mnt/home/huangbi4/PPI_project/new/results/PPI_data/STRING_v11_PPI_score.txt'

'''

from collections import defaultdict
import argparse


def Parser():
   parser=argparse.ArgumentParser()
   parser.add_argument('-i', help = 'the input file')
   parser.add_argument('-o', help = 'the output file')
   parser.add_argument('-s', help = 'the PPI score threshold', type = int, default = 100)
   parser.add_argument('-p', help = 'the ppi file')
   return parser.parse_args()


def StorePPI(ppi_file, score_threshold):
   ppi = defaultdict(int)   
   with open(ppi_file) as f:
      for line in f:
         tf1, tf2, score = line.strip().split('\t')
         if int(score) > score_threshold:
            ppi[(tf1,tf2)] = 1
            ppi[(tf2,tf1)] = 1
   return ppi


def FilterFile(input_file, output_file, ppi):
   f, w = open(input_file), open(output_file,'w')
   for line in f:
      line = line.strip().split('\t')
      tf_list = line[5:]
      store = []
      for tf_pair in tf_list:
         tf1, tf2 = tf_pair.split(',')
         if ppi[(tf1,tf2)]: store.append(tf_pair)
      if len(store) > 0: w.writelines('\t'.join(line[:5])+'\t'+'\t'.join(store)+'\n')
   f.close(), w.close()



def main():
   args = Parser()
   input_file, output_file, ppi_file, score = args.i, args.o, args.p, args.s
   ppi = StorePPI(ppi_file, score)
   FilterFile(input_file, output_file, ppi)
   return 0



###########
main()
