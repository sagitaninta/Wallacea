# Mutational Load
This section of the pipeline consists of four steps:

- Creating an input databese with conservation scores and priors
- Obtaining genotype likelihoods (using ANGSD[^1])
- Calculating posterior probabilities
- Estimating genetic load

Here we present our strategy to perform all calculation for a single individual but this can be easily parallelized to an arbitrary number of individuals in the dataset.
[^1]:Korneliussen, T.S., Albrechtsen, A. & Nielsen, R. ANGSD: Analysis of Next Generation Sequencing Data. BMC Bioinformatics 15, 356 (2014). https://doi.org/10.1186/s12859-014-0356-4 

## 1. Prior db


## 2. Genotype likelihood
This step requires three input, namely, a bam file (along with its full path), the sample ID and a list of genomic coordinate (one site per line).

```sh
mkdir genolik
chmod 770 scripts/angsd_genolik.sh
./scripts/angsd_genolik.sh path/to/bam sample_ID coordinate_file
```
Note that coordinates for ANGSD are 1-based, hence, these coordinates should match the `End` field in your `.bed` file (see `example_coord_file.txt` for correct formatting).

Assuming a tab separated file (`angsd_args_file.tsv`) containing the required arguments in the order shown above, this step can be parallelized to an arbitrary number of individuals as follows:

```sh
cat angsd_args_file.tsv | xargs -L1 -P0 ./scripts/angsd_genolik.sh
```

## 3. Posterior

## 4. Load

