This folder contains the ProTECT predicted enhancer-promoter interactions in GM12878 and K562 based on distance-aware q-values. The predictions are sorted based on the starting sites of the enhancers.
The files have 11 columns:
<br />1. Name of the chromosome.
<br />2. Start of the enhancer.
<br />3. End of the enhancer.
<br />4. Start of the promoter (the promoter is defined as the 2 kb region centered at the TSS).
<br />5. End of the promoter.
<br />6. The genomic distance between the center of the enhancer and the interacting TSS.
<br />7. The predicted posterior probability of being true enhancer-promoter interaction based on the random forest model.
<br />8. The p-value of the enhancer-promoter interaction prediction based on the permutation tests and the predicted posterior probability.
<br />9. The distance-aware q-value based on the pFDR procedure.
<br />10. The Ensembl gene id of the interacting gene.
<br />11. The name of the interacting gene.
