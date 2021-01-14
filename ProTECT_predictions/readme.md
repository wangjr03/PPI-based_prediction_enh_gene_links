This folder contains the ProTECT predicted enhancer-promoter interactions in GM12878 and K562 based on distance-aware q-values. The predictions are sorted based on the starting sites of the enhancers. 
The files have 11 columns:
	1. Name of the chromosome.
	2. Start of the enhancers.
	3. End of the enhancers.
	4. Start of the promoters (The promoter is defined as the 2 kb region centered at the TSS).
	5. End of the promoters.
	6. The Genomic distance between the center of the enhancers and the interacted TSS.
	7. The predicted posterior probability of being true enhancer-promoter interactions based on the random forest model.
	8. The p-value of enhancer-promoter interaction prediction based on the permutation tests and predicted posterior probability.
	9. The distance-aware q-values based on the pFDR procedure.
	10. The Ensembl gene id of the interacted gene.
	11. The name of the interacted gene.
