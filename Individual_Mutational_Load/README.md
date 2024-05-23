# Mutational Load
This section of the pipeline consists of three steps:

- Calculating genotype likelihoods (using ANGSD[^1])
- Obtaining posterior probabilities
- Calculating load

Here we present the strategy to complete each task for a single individual but this can be easily parallelized to an arbitrary number of individuals in the dataset.
[^1]:Korneliussen, T.S., Albrechtsen, A. & Nielsen, R. ANGSD: Analysis of Next Generation Sequencing Data. BMC Bioinformatics 15, 356 (2014). https://doi.org/10.1186/s12859-014-0356-4 

## 1. Genotype likelihood
This step requires three input, namely, a bam file (along with its full path), the sample ID and a list of genomic coordinate (one site per line).

```sh
mkdir genolik
chmod 770 angsd_genolik.sh
./angsd_genolik.sh path/to/bam sample_ID coordinate_file
```
