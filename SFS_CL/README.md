# Composite likelihood approach to differentiate site frequency spectrum
This section of the pipeline consists of two steps:

- Creating an input SFS using ANGSD[^1]
- Calculating composite likelihood from SFS based on SweepFinder[^2], then downsampling and bootstrapping

[^1]:Korneliussen, T.S., Albrechtsen, A. & Nielsen, R. ANGSD: Analysis of Next Generation Sequencing Data. BMC Bioinformatics 15, 356 (2014). https://doi.org/10.1186/s12859-014-0356-4 
[^2]:Nielsen, R., Williamson, S., Kim, Y., Hubisz, M.J., Clark, A.G., & Bustamante, C. Genomic scans for selective sweeps using SNP data. Genome Research 15, 1566–1575 (2005). https://doi.org/10.1101/gr.4252305.

## 1. Creating an input SFS

Within a set of predetermined set of sites within low, moderate, modifier, and high impact rating, we generate SFS from a set of BAM from each of the population following ANGSD pipeline.

First, we called the site allele frequency (SAF):
```
REF=<path to reference genome>
ANC=<path to ancestral fasta from set of outgroups mapped to the $REF>

sites=$1 # ANGSD sites file
pop=$2 # file containing list of BAM files representing one population
out=${sites#*107/} # getting the basename of the sites file without the path

angsd -nThreads 4 -bam $pop -anc $ANC -ref $REF -sites $sites -C 50 -minQ 30 -minMapQ 30 -dosaf 1 -GL 2 -out ${pop%.txt}_${out%.file}
```

Then we call the site frequency spectrum:
```
file=$1 # the saf.idx from the output of previous code

realSFS $file > ${file%.saf.idx}_unfolded.sfs
```

## 2. Composite likelihood

The script to calculate composite likelihood and conduct a bootstrapping of the composite likelihood is given in the `02_clr_bootstrap.R`. Scripts for the plots showing the normalised composite likelihood values are given in `03_plotSFS_CL_normAll.R`

