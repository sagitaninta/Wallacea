import pandas as pd
import numpy as np
import math
import sys
import os



# Defining input and output name
post_file=sys.argv[1]
sample=os.path.splitext(os.path.basename(post_file))[0]
outfile="load/"+sample+".load"

# Importing the set of load and posterior probability, one variant per row.
post_df = pd.read_csv(post_file, header=0, skip_blank_lines=True, sep="\t")

# Note that threshold for deleteriousness are hard coded. 
# please modify the following two lines accordingly
threshold1=1.5
threshold2=0.75

# This function calculates the load score for a single variant (one row in df)
def do_load(df,score,t):
	hom=df['Alt']+df['Alt']
	het="".join(sorted([df['Ref'],df['Alt']]))
	hom_load=df[score]*df[hom]
	tot_load=hom_load+df[score]*df[het]
	if t=='H':
		return(hom_load)
	else:
		return(tot_load)

# Applying the do_load function for different scores and load types
post_df['PP_hom']=post_df.apply(do_load, args=('phyloP','H'), axis=1)
post_df['PP_tot']=post_df.apply(do_load, args=('phyloP','T'), axis=1)
post_df['PC_hom']=post_df.apply(do_load, args=('phasCon','H'), axis=1)
post_df['PC_tot']=post_df.apply(do_load, args=('phasCon','T'), axis=1)

# Subsetting df for threshold
subset_PP=post_df[post_df['phyloP']>=threshold1]
load_pp_h=round(subset_PP['PP_hom'].sum()/len(subset_PP['PP_hom']),4)
load_pp_t=round(subset_PP['PP_tot'].sum()/len(subset_PP['PP_tot']),4)

# Subsetting df for threshold
subset_PC=post_df[post_df['phasCon']>=threshold2]
load_pc_h=round(subset_PC['PC_hom'].sum()/len(subset_PP['PC_hom']),4)
load_pc_t=round(subset_PC['PC_tot'].sum()/len(subset_PP['PC_tot']),4)

# Writing output file
with open(outfile, 'w') as of:
	header="phyloP_hom" + "\t" + "phyloP_tot" + "\t" + "phasCon_hom" + "\t" + "phasCon_tot" + "\n"
	of.write(header)
	res_line="\t".join([str(load_pp_h),str(load_pp_t),str(load_pc_h),str(load_pc_t)])+"\n"
	of.write(res_line)
