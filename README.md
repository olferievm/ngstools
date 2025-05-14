# ngstools
A collection of R scripts for RNAseq data analysis

Includes useful functions for fast access between project.

tpm - calculate TPM (Transcripts Per Million) from an edgeR object in R. based on the follow formula:

$`{TPM}_{ij} = \frac{\frac{\text{counts}_{ij}}{\text{length}_i}}{\sum_i \frac{\text{counts}_{ij}}{\text{length}_i}}`$

Where:
- i is gene 
- j is sample 
- $length_i$ is a gene length in kilobases (kb)

