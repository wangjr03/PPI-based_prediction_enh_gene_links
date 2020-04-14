Summary:
This package builds models based on the chromatin interaction data and ChIP-seq binding files, and predicts the probability of interaction between new input promoter-enhancer pairs.



Introduction:

This set of programs operate on chromatin interaction data (e.g. Hi-C, ChIA-PET) and ChIP-seq binding file. They can be used to
	1. Generate promoter-enhancer interaction pairs from the input chromatin interaction data
	2. Profile the transcription factors binding in the promoter-enhancer interaction pairs
	3. Train Random Forests model on the generated promoter-enhancer interaction pairs, and predict the probability of interaction for a new chromatin region pairs


Dependencies:
   1. Python3
   2. Python packages: scipy, sklearn, matplotlib, numpy, pandas


Files:
   1. generate_promoter_list.py: generate cell-type specific promoter list
   2. get_enhancer_promoter_pair.py: generate promoter-enhancer interaction pairs from the input chromatin interaction data
   3. generate_negative_group.py: generate false promoter-enhancer interaction pairs (negative set). This negative set will be used to train the Random Forests model
   4. TF_overlapping.py: profile the transcription factors binding in the promoter-enhancer interaction pairs
   5. filter_based_on_ppi_score_remove_none_TFpairs.py: remove the transcription factor pairs in the promoter-enhancer pairs that are not supported by protein-protein interaction data
   6. random_forest.py: train a Random Forests model using positive set and negative set, and do predictions on new input data




General Usage:
   1: generate_promoter_list.py
           -c: the cell type name
           -f: define the promoter flanking region, +/- flanking_region is considered as the promoter region
           -o: the output promoter list file
           -cell_type: the cell type list file, default = Cell_types_RNAseq.bed
           -expression: the gene expression file, default = RPKM_all_gene_56epigenomes_simple_select_2.bed 
           -gene_name: the gene name file, default = RPKM_all_gene_name_select_2.bed
           -gene_pos: the gene position file, default = Gene_annotation.bed

   2: get_enhancer_promoter_pair.py
           -i: the chromatin interaction data
           -o: the output enhancer-promoter pairs

   3: generate_negative_group.py
           -i: the true enhancer-promoter pairs (positive set).
           -o: output the generated false promoter-enhancer interaction pairs (negative set). The distance distributions of the positive set and the negative set are the same
           -n: downsample n enhancers from the whole enhancer list set
           -b: the bin size to calculate the distance distribution (unit:bp)
           -t: the relative size of negative set to the positive set

   4: TF_overlapping.py
           -i: the input enhancer-promoter pairs (positive set or negative set)
           -o: the output TF profiles in enhancer-promoter pairs
           -p: the path to the transcription factor peak files
           -c: the column of chromosome, enhancer start site, enhancer end site, promoter start site, promoter end site

   5: filter_based_on_ppi_score_remove_none_TFpairs.py
           -i: the input enhancer-promoter pairs with transcription factor pairs features
           -o: output enhancer-promoter pairs with transcription factor pairs features supported by PPI data
           -s: the PPI score threshold
           -p: the PPI data

   6: random_forest.py
           -i: input the positive set and the negative set
           -l: the corresponding label of the input sets (1: positive set, -1: negative set)
           -k: k-fold Cross-Validation
           -s: the features in the Random Forests model and the corresponding feature importances in the k-fold Cross-Validation
           -activity: the cell line. This parameter is used when including the gene activity, enhancer activity, enhancer-gene correlation features
           -d: if to include the distance features in the Random Forests model
           -n: the number of decision trees in the Random Forests model 
           -newin: the new input region pair data
           -newout: the prediction on whether there are interactions between the input region pairs
           -o: the saved ROC curve image
