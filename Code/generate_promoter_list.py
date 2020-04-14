'''
-c: the cell line key words

-o: the output file name

-f: the flanking size of TSS to get a promoter, promoter = TSS +/- flanking_size, default = 1000, type = int

-cell_type: the cell type file, default = '/Users/hbb/Dropbox/PPI_project/Data/gene_activity/Cell_types_RNAseq.bed'

-expression: the gene expression file, default = '/Users/hbb/Dropbox/PPI_project/Data/gene_activity/RPKM_all_gene_56epigenomes_simple_select_2.bed'

-gene_name: the gene name file, default = '/Users/hbb/Dropbox/PPI_project/Data/gene_activity/RPKM_all_gene_name_select_2.bed'

-gene_pos: the gene position file, default = '/Users/hbb/Dropbox/PPI_project/Data/gene_activity/Gene_annotation.bed'
'''


import numpy as np
import re
import pandas as pd
import argparse

def Parser():
   parser=argparse.ArgumentParser('')
   parser.add_argument('-c', help = 'the cell line key words')
   parser.add_argument('-o', help = 'the output file name')
   parser.add_argument('-f', help = 'the flanking size of TSS to get a promoter, promoter = TSS +/- flanking_size', default = 1000, type = int)
   parser.add_argument('-cell_type', help = 'the cell type file', default = '/Users/hbb/Dropbox/PPI_project/Data/gene_activity/Cell_types_RNAseq.bed')
   parser.add_argument('-expression', help = 'the gene expression file', default = '/Users/hbb/Dropbox/PPI_project/Data/gene_activity/RPKM_all_gene_56epigenomes_simple_select_2.bed')
   parser.add_argument('-gene_name', help = 'the gene name file', default = '/Users/hbb/Dropbox/PPI_project/Data/gene_activity/RPKM_all_gene_name_select_2.bed')
   parser.add_argument('-gene_pos', help = 'the gene position file', default = '/Users/hbb/Dropbox/PPI_project/Data/gene_activity/Gene_annotation.bed')
   return parser.parse_args()




def StoreCellLine(cell_type_file):
   return pd.read_csv(cell_type_file, sep = '\t')['Universal_Human_Reference'].tolist()



def CellLineKeyWord(cell_type_list, input_cell_line):
   for i in cell_type_list:
      if re.search(input_cell_line.upper(), i.upper()): return i
   return ''


def StoreExpression(gene_expression_file, cell_type_list, target_cell_type):
   return pd.read_csv(gene_expression_file, names = cell_type_list, sep = '\t')[target_cell_type].tolist()



def StoreGeneList(gene_name_file):
   return pd.read_csv(gene_name_file, sep = '\t', header = None)[0].tolist()


def FilterGene(gene_list, expression_profile):
   return [gene_list[i] for i in range(len(gene_list)) if expression_profile[i]]

def GetPromoter(gene_pos_file, active_gene_list, flanking_size):
   store = set()
   with open(gene_pos_file) as f:
      for line in f:
         ensembl_id, chrom, start, end, strand, type, gene_symbol = line.strip().split('\t')[:7]
         if ensembl_id in active_gene_list:
            tss = int(start) if strand == '1' else int(end)
            store.add(('chr'+chrom, tss-flanking_size, tss+flanking_size, gene_symbol))
   return store



def Output(promoter_list, output_file):
   with open(output_file,'w') as w:
      for a_promoter in promoter_list: w.writelines('\t'.join(np.array(a_promoter).astype(str))+'\n')
   return 0



def main():
   args = Parser()
   cell_type_file, cell_type_key, expression_file, gene_name_file, gene_pos_file, flanking_size, output_file = args.cell_type, args.c, args.expression, args.gene_name, args.gene_pos, args.f, args.o 
   cell_line_list = StoreCellLine(cell_type_file)
   cell_line = CellLineKeyWord(cell_line_list, cell_type_key)
   if cell_line == '': return 'Error, cell type not found'
   expression_profile = StoreExpression(expression_file, cell_line_list, cell_line)
   gene_list = StoreGeneList(gene_name_file)
   active_gene_list = FilterGene(gene_list, expression_profile)
   promoter_list = GetPromoter(gene_pos_file, active_gene_list, flanking_size)   
   Output(promoter_list, output_file)
   return 0 







##########
main()
