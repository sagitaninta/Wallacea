import pandas as pd
import numpy as np
import math
import sys
import os



# Defining input and output name
post_file=sys.argv[1]
sample=os.path.splitext(os.path.basename(post_file))[0]
outfile="sift/load/"+sample+".load"

# Importing the set of load and posterior probability, one variant per row.
post_df = pd.read_csv(post_file, header=0, skip_blank_lines=True, sep="\t")
threshold=0.05
# this function calculates the load score for a single variant (one row in df)
def do_load(df,score,t):
	hom=df['Alt']+df['Alt']
	het="".join(sorted([df['Ref'],df['Alt']]))
	hom_load=(1-df[score])*df[hom]
	tot_load=hom_load+(1-df[score])*df[het]
	if t=='H':
		return(hom_load)
	else:
		return(tot_load)

# applying the do_load function for different scores and load types
post_df['sift_hom']=post_df.apply(do_load, args=('sift','H'), axis=1)
post_df['sift_tot']=post_df.apply(do_load, args=('sift','T'), axis=1)

# subsetting df for threshold
subset_sift=post_df[post_df['sift']<=threshold]
load_h=round(subset_sift['sift_hom'].sum()/len(subset_sift['sift_hom']),4)
load_t=round(subset_sift['sift_tot'].sum()/len(subset_sift['sift_tot']),4)

with open(outfile, 'w') as of:
	header="sift_hom" + "\t" + "sift_tot" + "\n"
	of.write(header)
	res_line="\t".join([str(load_h),str(load_t)])+"\n"
	of.write(res_line)
