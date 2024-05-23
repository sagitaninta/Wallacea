import pandas as pd
import itertools
import sys

# These two functions creates a prior probability mass function based on
# the alleles segregating at each site and how many they are
# They relay on the variable genotypes defined below.
# That variable serves as keys of the dictionary {Genotype: Probability}

genotypes=list(itertools.combinations_with_replacement(["A","C","G","T"],2))
genotypes.sort()

def alleles(a,b):
	allele_list=list(a)+b.split(",")
	allele_list.sort()
	genotype_keys=list(itertools.combinations_with_replacement(allele_list,2))
	return(genotype_keys)

def prior(keys):
	prior_dict=dict.fromkeys(genotypes,0.01)
	if len(keys)==3:
		seg_dict=dict.fromkeys(keys,0.31)
	elif len(keys)==6:
		seg_dict=dict.fromkeys(keys,0.16)
	else:
		seg_dict=dict.fromkeys(keys,0.1)
	prior_dict.update(seg_dict)
	prior_values=[prior_dict[key] for key in prior_dict.keys()]
	return(prior_values)


# Importing simplified vcf containing one variant per row
snp_data_file = sys.argv[1]
outfile = sys.argv[2]
snp_df = pd.read_csv(snp_data_file, header=None, skip_blank_lines=True, sep="\t")
snp_df.columns = ["Chromosome","start", "end", "Variant_ID", "Ref", "Alt", "phasCon","phyloP"]

# Creating prior data frame and merging it with relevant column of the snp data frame
prior_df=snp_df.apply(lambda row: prior(alleles(row.Ref, row.Alt)), axis=1, result_type="expand")
prior_df.columns = ["AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"]
variants=snp_df[["Variant_ID","Ref", "Alt", "phasCon", "phyloP"]].join(prior_df)
variants.to_csv(outfile, header=True, index=False, sep="\t")
