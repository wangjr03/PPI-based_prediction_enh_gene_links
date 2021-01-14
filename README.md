# ProTECT
Predict long-range enhancer-promoter intearctions based on PPIs.
## Summary
We developed ProTECT (i.e. PROtein-protein interactions of Transcription factors predicting Enhancer Contacts with Target genes) to predict long-range enhancer-promoter interactions using TF PPI as features. A robust random forest model is trained based on significant experimental chromatin interactions,i.e. Hi-C, and applied to the whole genome-wide to make highly confident predictions. <\p>

## Introduction
The ProTECT algorithm takes multi-omics data as inputs, including the enhancer activities, gene activities, TF Chip-seq narrow peaks and PPIs, to train a random forest model. Novel community detection-based feature dimension reduction and feature selection methods are applied to reduce the number of features and improve the robustness and accuracy of the ProTECT. In the meanwhile, we evaluate the accuracy of the ProTECT with a rigorous genomic-bin split cross-validation method to control the confounding factors. By applying the ProTECT algorithm on the potential enhancer-promoter interactions in the whole genome-wide, a set of highly confident enhancer-promoter interactions is prioritized and used for the downstream analysis.


## Dependencies
The implementation of the ProTECT is based on `Python 3.6` and `R 3.5.1`. It depends on 6 Python packages (`scipy`, `scikit-learn`, `numpy`, `pandas`,`pickle` and `re`) and 2 R packages (`igraph` and `expm`).

## Input data: resources and formats
The ProTECT takes 8 sets of data as inputs: (1) significant experimental chromatin interactions, (2) contact domain annotations, (3) enhancer annotations, (4) gene annotations, (5) the enhancer activity profiles, (6) the gene activity profiles, (7) TF ChIP-seq peaks and (8) protein-protein interaction (PPI) datasets. For the convenience of the user,  enhancer annotations, gene annotations, the enhancer activity profiles and the gene activity profiles are pre-calculated and provided. Descriptions of the three required user-provided data are listed below:
Significant experimental chromatin interactions: The significant experimental chromatin interactions can be provided by Hi-C , ChIA-PET  and Capture C . The chromatin interactions should be in tab separated file with five columns:
	| col | abbrv. | type | description |
	| --- | --- | --- | --- |
	| 1 | chr | string | Name of the chromosome |
	| 2 | frag1.start | int | Fragment 1 start |
	| 3 | frag1.enh | int | Fragment 1 enh |
	| 4 | frag2.start | int | Fragment 2 start |
	| 5 | frag2.end | int |Fragment 2 end |
2. Contact domain annotation: The contact domains represent densely self-interacting genome regions. The contact domain can be detected by applying computational models, e.g. Arrowhead , on chromatin contact maps. For example, contact domains based on Hi-C contact maps can be downloaded from GEO with GSE63525 . The contact domain annotations should be tab-separated, with three columns:
	| col | abbrv. | type | description |
	| --- | --- | --- | --- |
	| 1 | chr | string | Name of the chromosome |
	| 2 | domain.start | int | Contact domain start |
	| 3 | domain.end | int | Contact domain end |

3. TF ChIP-seq peaks: The narrow peak files of TF ChIP-seq can be downloaded from the ENCODE consortia. As the quality control, we applied three criterions to filter the TF ChIP-seq narrow peak files. The TF ChIP-seq narrow peak files with the best quality are selected for each TF using the following three criteria:
	a. TF ChIP-seq peak files for treated transcription factors are removed.
	b. TF ChIP-seq peak files generated by paired-end experiments are preferred if available.
	c. FRiP (Fraction of Reads in Peaks) score is calculated for each TF ChIP-seq replicate. The TF ChIP-seq peak file with the highest averaged FRiP score is selected.
The format of the TF ChIP-seq narrow peak files follows the standard definition of the narrow peak file.

## Pre-calculated data
For users’ convenience, four sets of data have been pre-calculated. The user could also use their own datasets by replacing those files.
1. Gene annotations: The gene annotation with GENCODE V17 has been integrated with the program. The promoter is defined as the +/- 1kb around the transcriptional start sites (TSS). The Gene annotations should be in the following format.
	| col | abbrv. | type | description |
	| --- | --- | --- | --- |
	| 1 | ensg.id | char | Ensembl gene id |
	| 2 | chr | char | Name of the chromosome |
	| 3 | gene.start | int | Gene body start |
	| 4 | gene.end | int | Gene body end |
	| 5 | strand | int | Strand of the gene, 1 for ‘+’ and -1 for ‘-’ |
	| 6 | gene.type | int | Type of the gene, i.e. protein_coding |
	| 7 | gene.name | char | Name of the gene, i.e. GATA1|
	| 8 | HGNC description | string | Description of the gene based on HGNC annotation |
	
	
2. Enhancer activity matrix: The enhancer activities are quantified by cell-type specific genome-wide coverage epigenomic datasets, i.e. histone modifications, DNase-seq and ATAC-seq. The cell-type specific enhancer activities are summarized into a matrix, where rows present enhancers, and columns represent cell-types.
3. Gene expression matrix: The gene expressions are quantified by RNA-seq data, i.e. RPKM. The cell-type specific gene expressions are summarized into a matrix, where rows present genes, and columns represent cell-types.
4. Protein-protein interaction (PPI) datasets: PPI datasets can be downloaded from STRING database. To remove the low-quality PPIs, we only used the PPIs with the confidence score greater than 0 in the ‘Experiments’ category. A filtered PPI dataset based on STRING V11 has been integrated into the program.

## Description of scripts: command lines
The ProTECT software consists of 6 sequential scripts. A detailed description of each piece is provided.
1. Training_sample_generation.R: This step is used to generate positive training sets and a balanced negative training set with multiple confounding factors controlled. 
	**Inputs**: It takes significant experimental chromatin interactions, contact domain annotations as inputs. 
	**Outputs**: a list of enhancer-promoter interactions and their labels, i.e. 1 for positive sets and -1 for negative sets.
	**Command line usage**: Rscript Training_sample_generation.R `<path to significant chromatin interactions>` `<path to contact domain annotations>`
2. Generate_ProTECT_feature.py: This step is used to generate features used by ProTECT, including enhancer activities, gene expressions, enhancer-promoter activity correlations, genomic distances and PPI features.
	**Inputs**: It takes TF ChIP-seq narrow peak files as inputs.
	**Outputs**: a feature matrix, where each row is one enhancer-promoter interaction and each column is one feature.
	**Command line usage**: python Generate_ProTECT_feature.py -i `<path to the training samples>` -c `<column index of the cell line in the gene/enhancer activity matrix, i.e. 53 for the GM12878 cell line>` -d `(indicates the genomic distance is reported)` -s `<path to the PPI file>`

3. Discover_PPI_module.R: This step is used to detect a two -layer hierarchical PPI networks.
	**Inputs**: It takes the PPI data as inputs.
	**Outputs**: Membership of the TFs to the hierarchical PPI module.
	| col | abbrv. | type | description |
	| --- | --- | --- | --- |
	| 1 | TF_name | string | Name of the TF |
	| 2 | S_module.index | int | Module index for S-module |
	| 3 | L-module.index | int | Module index for L-module |	

	**Command line usage**: Rscript Discover_PPI_module.R `<path to the PPI data>` `<threshold of the confidence score>`
	
4. ProTECT.py: This script is used to do feature dimension reduction and train a random forest model based on the generated feature.
	**Inputs**: It takes the outputs of step 1,2,3 as the inputs.
	**Outputs**: a trained model stored in pickle format and a text summary of cross-validation results.
	**Command line usage**: python ProTECT.py -o `<output directory>`

5. For genome-wide application, users should provide the potential enhancer-promoter interactions. A recommended method is to use the bedtools window function. Following is an example: bedtools window -a <promoter annotation> -b <enhancer annotation> -w <maximum distance between the enhancer and the promoter> > <outputpath of potential enhancer-promoter interactions>. To generate the feature matrix based on individual TF-TF pairs on these potential pairs, we can repeat step 2 and the output is the desired feature matrix (raw feature matrix). 

6. Reformat_feature_predict.py: This script uses the trained random forest model to predict significant enhancer-promoter interactions from the whole pool. This script has two major steps. The first step aims to reformat the TF-level PPI feature matrix generated in step 5 into the module-level PPI features defined by step 4. The second step takes the reformatted feature matrix as input and assigns a probability to each enhancer-promoter interaction using the training random forest model.
	**Inputs**: the TF-level PPI feature matrix generated in step 5.
	**Outputs**: the reformatted feature matrix defined by feature engineering procedures in step 4 and a file containing the predictive probability for each potential enhancer-promoter interaction.
	**Command line usage**: python Reformat_feature.py -i `<path to the raw feature matrix>`

7. Reformat_feature_shuffle.py: This script is used to generate the null probability distribution for potential enhancer-promoter interactions. In this step, each feature is permuted randomly across all samples. The permuted features are then used as the input of the trained random forest model to generate the null probabilities. The input and usage are the same as step 6. 

8. pFDR_procedure.R: This script uses the pFDR procedure to integrate the genomic distance information. The first step is to calculate a p-value for each enhancer-promoter interaction based on the predictive probability from step 6 and null predictive probabilities step 7. The second step is to calculate a q-value using the pFDR procedure based on the genomic distances and p-values.
	**Inputs**: the result of step 6 and step 7.
	**Outputs**: A file containing q-values and statistical significant enhancer-promoter interactions given the q-value threshold.
	**Command line usage**: Rscript pFDR_procedure.R `<path to the predictive probability>` `<path to the null predictive probability> <q-value threshold>`
	
	




	
