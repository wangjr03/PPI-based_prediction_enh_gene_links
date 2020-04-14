import argparse
import numpy as np
import random
import time
from collections import defaultdict

parser=argparse.ArgumentParser('')
parser.add_argument('-i',help='input file',required=True)
parser.add_argument('-o',help='output file',required=True)
parser.add_argument('-n', help = 'the number of enhancers downsampled from enhancer_coords_select.bed file', type = int, default = 60000)
parser.add_argument('-b', help = 'bin size for distance distribution', type = int, default = 10000)
parser.add_argument('-t', help = 'negative set is t times of the positive set', type = int, default = 1)
args=parser.parse_args()


random.seed(int(time.time()))


########first get the frequency of every group
##1. store the distance information for positive group
enhancer, promoter, positive_group = defaultdict(list), defaultdict(list), defaultdict(list)
positive_dist, positive_pool = [], []


input_file = args.i
with open(input_file) as f:
   for line in f:
      line=np.array(line.strip().split('\t'))
      chrom=line[0]
      enhancer_start, enhancer_stop, promoter_start, promoter_stop = line[[1,2,4,5]].astype(int)
      enhancer_middle=enhancer_start+(enhancer_stop-enhancer_start)/2
      promoter_middle=promoter_start+(promoter_stop-promoter_start)/2
      d=int(abs(promoter_middle-enhancer_middle))
      enhancer[chrom].append((enhancer_start,enhancer_stop))
      promoter[chrom].append((promoter_start,promoter_stop))
      positive_group[chrom].append((enhancer_start,enhancer_stop,promoter_start,promoter_stop))
      positive_dist.append(d)

print('finish store positive set\n')


downsample_enhancer_n = args.n
index = sorted(random.sample(range(697876),downsample_enhancer_n))
count, point = 0, 0
with open('/Users/hbb/Dropbox/PPI_project/Data/enhancer/enhancer_coords_select.bed') as f:
   for line in f:
      if point == downsample_enhancer_n: break
      if count == index[point]:
         chrom, start, end =line.strip().split('\t')[:3]
         enhancer[chrom].append((int(start),int(end)))
         point +=1
      count +=1


print('finish downsample enhancer positions\n')



##2. generate negative group pools
negative_pool=[]
for chrom in enhancer:
   for enhancer_start, enhancer_stop in enhancer[chrom]:
      for promoter_start, promoter_stop in promoter[chrom]:
         link=(enhancer_start,enhancer_stop,promoter_start,promoter_stop)
         if (link in positive_group[chrom]): pass
         else: negative_pool.append((chrom,enhancer_start,enhancer_stop,promoter_start,promoter_stop))


print('finish generate negative pool\n')




##3. get the min, max, gap information, calculate the counts in every region for positive groups
bin_size = args.b
min_val = int(min(positive_dist)/bin_size) * bin_size
max_val = int(max(positive_dist)/bin_size+1) * bin_size

count_array, range_array = np.histogram(positive_dist, bins = range(min_val, max_val+bin_size, bin_size))
positive_count = dict(zip(range_array, count_array))

print('finish the positive distance distribution\n')

 
##4. distribute the negative groups according to its distance
negative_distribution=defaultdict(list)
for chrom, enhancer_start, enhancer_stop, promoter_start, promoter_stop in negative_pool:
   enhancer_middle=enhancer_start+(enhancer_stop-enhancer_start)/2
   promoter_middle=promoter_start+(promoter_stop-promoter_start)/2
   d=int(abs(promoter_middle-enhancer_middle))
   if d < min_val or d > max_val: pass
   else:
      key=int((d-min_val)/bin_size)*bin_size + min_val
      negative_distribution[key].append((chrom,enhancer_start,enhancer_stop,promoter_start,promoter_stop))

print('finish assign each negative link to its corresponding list\n')


##6. sampling from the positive and negative pools
w= open(args.o,'w')
times = args.t
for key in positive_count.keys():
   needed_count=positive_count[key]*times
   pool_count=len(negative_distribution[key])
   if needed_count == 0: pass
   elif pool_count == 0: print('No negative pairs in '+str(key)+'-'+str(key+bin_size)+'\n')
   elif  needed_count >= pool_count:
      print('No enough negative pairs in '+str(key)+'-'+str(key+bin_size)+'\n')
      for a_pair in negative_distribution[key]: w.writelines('\t'.join(np.array(a_pair).astype(str))+'\n')
   else:
      downsample = random.sample(negative_distribution[key],needed_count)
      for a_pair in downsample: w.writelines('\t'.join(np.array(a_pair).astype(str))+'\n')

print('finish output the negative set\n')

