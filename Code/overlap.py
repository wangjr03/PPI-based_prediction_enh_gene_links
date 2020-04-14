from collections import defaultdict
class overlap():
   def modify_dict(self, input_dict):
      store = {}
      for a_key in input_dict: store[a_key] = list(input_dict[a_key])
      return store
   
   def do_overlap(self,input1,input2): 
      overlapping = {}
      chromosome = set(input1.keys()) & set(input2.keys())
      for chrom in chromosome:
         overlapping[chrom] = defaultdict(set)
         input1_peak, input2_peak = sorted(list(input1[chrom])), sorted(list(input2[chrom]))
         pointer, total_peak = 0, len(input2_peak)
         for start1, end1 in input1_peak:
            if end1 < input2_peak[pointer][0]: pass
            elif start1 > input2_peak[pointer][1]:
               if pointer == total_peak-1: pass
               else:
                  while (start1 > input2_peak[pointer][1] and pointer < total_peak-1): pointer+=1
                  if end1 < input2_peak[pointer][0]: pass
                  elif start1 > input2_peak[pointer][1]: pass
                  else:
                     start2, end2 = input2_peak[pointer]
                     overlapping[chrom][(start1,end1)].add((start2,end2))
                     pointer1 = pointer
                     while (pointer1 < total_peak-1  and  end1 > input2_peak[pointer1+1][0]):
                        pointer1 += 1
                        start2, end2 = input2_peak[pointer1]
                        overlapping[chrom][(start1,end1)].add((start2,end2))
            else:
               start2, end2 = input2_peak[pointer]
               overlapping[chrom][(start1,end1)].add((start2,end2))
               pointer1 = pointer
               while (pointer1 < total_peak-1  and  end1 > input2_peak[pointer1+1][0]):
                  pointer1 += 1
                  start2, end2 = input2_peak[pointer1]
                  overlapping[chrom][(start1,end1)].add((start2,end2))
      return overlapping


         
