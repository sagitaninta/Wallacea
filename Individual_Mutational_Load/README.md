# Mutational Load
This simple python based pipeline consists of four steps:

- Creating an input database with conservation scores and priors
- Obtaining genotype likelihoods (using ANGSD[^1])
- Calculating posterior probabilities
- Estimating genetic load

The first three steps essentially apply the Baesyan framework described in Plassais et al. (2022)[^2] and extend it to a small panel of SNPs. Here we present our strategy to perform all calculation for a single individual but this can be easily parallelized to an arbitrary number of individuals in the dataset. The pipeline has been designed for `python 3.11` and requires the following python libraries: `pandas`, `itertools`, `numpy`, along with `math`, `sys`, and `os`. For simplicity, we recommend to install these in a separate conda envirnment to avoid conflict and dependency issues.


[^1]:Korneliussen, T.S., Albrechtsen, A. & Nielsen, R. (2014). ANGSD: Analysis of Next Generation Sequencing Data. BMC Bioinformatics 15:356. https://doi.org/10.1186/s12859-014-0356-4 
[^2]: Plassais, J., vonHoldt, B.M, Parker, H.G., Carmagnini, A., et al. (2022). Natural and human-driven selection of a single non-coding body size variant in ancient and modern canids. Current Biology, 32:889-897.

## 1. Prior db
The pipeline requires a .bed file (Browser Extensible Data format) as input with the following fields: 

- Chromosome
- Starting position
- Ending position
- Variant ID
- Reference allele
- Alternative allele
- phasCons score
- phyloP score

Note that coordinate are zero-based (see `example-input_db.bed` for correct formatting).

The script `prior_creator.py` will parse the input database and output a .tsv file containing the variant ID, the allele segregating at each site and a list of informed prior probability values for each possible genotype (see `example_prior_db.tsv`). The script requires `python 3.9` and  
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

