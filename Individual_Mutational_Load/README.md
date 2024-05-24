# Mutational Load
This simple python based pipeline consists of four steps:

- Creating an input database with conservation scores and priors
- Obtaining genotype likelihoods (using ANGSD[^1])
- Calculating posterior probabilities
- Estimating genetic load

The first three steps essentially apply the Baesyan framework described in Plassais et al. (2022)[^2] to a small panel of SNPs. 

Here we present our strategy to perform all calculation for a single individual but this can be easily parallelized to an arbitrary number of individuals in the dataset. 

The pipeline has been designed for `python 3.11` and requires the following python libraries: `pandas`, `itertools`, `numpy`, along with `math`, `sys`, and `os`. For simplicity, we recommend to install these in a `conda` envirnment (hereafter termed `load`) to avoid conflict and dependency issues.


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

Note that coordinate are zero-based (see `example_input_db.bed` for correct formatting).

The script `prior_creator.py` will parse the input database and output a .tsv file containing the variant ID, the allele segregating at each site and a list of informed prior probability values for each possible genotype (see `example_priors_scores_db.tsv`). 

```sh
conda activate load
python input_db.bed priors_scores_db.tsv
```

## 2. Genotype likelihood
This step requires three input, namely, a bam file, the sample ID and a file containing the list of genomic coordinate (one site per line).

```sh
mkdir genolik
chmod 770 scripts/angsd_genolik.sh
./scripts/angsd_genolik.sh path/to/bam sample_ID coordinate_file
```
Note that coordinates for ANGSD are 1-based, hence, these coordinates should match the `Ending position` (3<sup>rd</sup> field) in your `input_db.bed` file (see `example_coord_file.txt` for correct formatting).

Assuming a tab separated file (`angsd_args_file.tsv`) containing the required arguments in the order shown above, this step can be parallelized to an arbitrary number of individuals as follows:

```sh
cat angsd_args_file.tsv | xargs -L1 -P0 ./scripts/angsd_genolik.sh
```
This will generate two files for each individuals in the `genolik` directory.

## 3. Posterior
After completing the previous steps we can run:
```sh
mkdir post
conda activate load
python posterior_estimator.py priors_scores_db.tsv ind_genotype_likelihood.glf.gz ./post/ind_posterior_file.post
```
The `posterior_estimator.py` script will calculate posterior probability for each site in the `output_db.tsv` for a given individual. 

As before this step can be easily parallelized by running:
```sh
mkdir post
ls ./genolik/*.glf.gz > posterior_list.txt
chmod 770 do_post.sh
conda activate load
cat posterior_list.txt | xargs -L1 -P0 ./do_post.sh
```
Where `do_post.sh` is simply a bash wrapper for the `posterior estimator.py` script.

## 4. Load
Finally to calculate mutational load we can run:
```sh
python load_estimator.py ./post/ind_file.post
```
which can be parallelized following the same strategy described above and using the `load.sh` wrapper.

## Calculating genetic load using SIFT scores
For legacy reasons we have created a separate pipeline for calculating genetic load using SIFT scores. This consists of the same four steps detailed above and employs variants of the python and bash scripts mentioned in the previous section with the prefix `sift`.
