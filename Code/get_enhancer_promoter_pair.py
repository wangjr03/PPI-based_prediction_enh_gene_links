'''
-i: the input file

-o: the output file

Output format:
enhancer_chrom	enhancer_start	enhancer_end	promoter_chrom	promoter_start	promoter_end
sep = '\t'

'''

import re
import argparse
import numpy as np


def Parser():
   parser=argparse.ArgumentParser('')
   parser.add_argument('-i', help = 'the input file')
   parser.add_argument('-o', help = 'the output file')
   return parser.parse_args()

  
def GetPair(input_file):
   pairs = []
   with open(input_file) as f:
      for line in f:
         if re.search('enhancer',line) and re.search('promoter',line):
            enhancer, promoter = EnhancerPromoterList(line)
            pairs.extend(EnhancerPromoterPair(enhancer, promoter))
   return set(pairs) 


def EnhancerPromoterList(line):
   enhancer, promoter = [], []
   line = line.strip().split('\t')
   for i in line[6:]:
      if re.search('promoter',i): 
         chrom, start, end = re.split('[(),]',i)[2:5]
         promoter.append((chrom, int(start), int(end)))
      else:
         chrom, start, end = re.split('[(),]',i)[2:5]
         enhancer.append((chrom, int(start), int(end)))
   return enhancer, promoter


def EnhancerPromoterPair(enhancer, promoter):
   pairs = []
   for a_enhancer in enhancer:
      for a_promoter in promoter:
         if not CheckSpaceOverlapping(a_enhancer, a_promoter):
            pairs.append((a_enhancer[0],a_enhancer[1],a_enhancer[2],a_promoter[0],a_promoter[1],a_promoter[2]))
   return pairs


def CheckSpaceOverlapping(region1, region2):
   chrom1, start1, end1 = region1
   chrom2, start2, end2 = region2
   if end2 < start1 or start2 > end1: return 0 ###not overlapping
   return 1



def Output(input_pair, output_file):
   with open(output_file,'w') as w:
      for a_pair in input_pair: w.writelines('\t'.join(np.array(a_pair).astype(str))+'\n')
   return 0





def main():
   args = Parser()
   pairs = GetPair(args.i)
   Output(pairs, args.o)
   return 0





########
main() 
