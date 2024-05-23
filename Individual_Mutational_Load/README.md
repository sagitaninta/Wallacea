# This section of the pipeline consists of three steps:

- Calculating genotype likelihoods (using ANGSD)
- Obtaining posterior probabilities
- Calculating load

Here we present the strategy to complete each task for a single individual but this can be easily parallelized to an arbitrary number of individuals in the dataset.

# 1. Genotype likelihood
This step requires three input, namely, a bam file (along with its full path), the sample ID and a list of genomic coordinate (one site per line).
