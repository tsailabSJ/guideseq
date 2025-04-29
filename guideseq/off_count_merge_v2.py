#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import glob
from joblib import Parallel, delayed
import re
import regex
import sys
import gc
import swifter
import warnings
import os
warnings.filterwarnings("ignore")

def merge_sites(d):
    d=d.sort_values(['distance_between_cas9_cut_and_most_frequent_guideseq_read_start'])
    dp = d.head(n=1)
    for c in [ 'total_guideseq_reads', 'total_control_reads', 'total_guideseq_non_dsODN_reads', 'total_control_non_dsODN_reads']:
        dp.loc[:, c] = d[c].sum()
    for c in [ 'nuclease.rev.plus', 'nuclease.fwd.plus', 'nuclease.fwd.minus', 'nuclease.rev.minus', 'control.rev.plus', 'control.fwd.plus', 'control.fwd.minus', 'control.rev.minus', 'nuclease.mispriming.plus', 'nuclease.mispriming.minus', 'control.mispriming.plus', 'control.mispriming.minus']:
        dp.loc[:,f"total.{c}"] = d[c].sum()
    dp.loc[:,f'merged_start_pos'] = ", ".join(d['end'].astype(str).tolist())
    dp.loc[:,"Number_assignments"] = d['Number_assignments'].max()
    return dp
# label = sys.argv[1]
inFile = sys.argv[1]
ncore = int(sys.argv[2])
df = pd.read_csv(inFile,sep="\t")
group_sizes = df.groupby(['#chr', 'start', 'end']).size()

# Map the group sizes back to the original DataFrame
df['Number_assignments'] = df.set_index(['#chr', 'start', 'end']).index.map(group_sizes)



df_list = [group for _, group in df.groupby(['#chr','start_absolute','end_absolute','strand','on_target_sequence'])]

out_list = Parallel(n_jobs=ncore,verbose=10,batch_size=10)(delayed(merge_sites)(d) for d in df_list)
df = pd.concat(out_list)
del out_list
gc.collect()

label = inFile.replace(".identified.raw.tsv","")
df.index = df['#chr']+":"+df.start_absolute.astype(str)+"-"+df.end_absolute.astype(str)+"_"+df.on_target_sequence

# slow

# identify same on-target self-overlap and merge
tmp_file=f"{label}.tmp.bed"
tmp_out = f"{label}.tmp.out"
df[['#chr','start_absolute','end_absolute','on_target_sequence']].to_csv(tmp_file,sep="\t",header=False,index=False)
command = f"bedtools intersect -a {tmp_file} -b {tmp_file} -wo > {tmp_out}"
os.system(command)
tmp = pd.read_csv(tmp_out,sep="\t",header=None)
tmp=tmp[tmp[8]>1] # only consider overlap>1bp 
tmp=tmp[tmp[3]==tmp[7]] # only consider same on-target
print (tmp)
tmp['S1'] = tmp[0]+":"+tmp[1].astype(str)+"-"+tmp[2].astype(str)+"_"+tmp[3]
tmp['S2'] = tmp[4]+":"+tmp[5].astype(str)+"-"+tmp[6].astype(str)+"_"+tmp[7]
tmp = tmp[tmp['S1']!=tmp['S2']]
self_overlap_dict={}
for i,r in tmp.iterrows():
	if r.S1 in self_overlap_dict:
		self_overlap_dict[r.S1].append(r.S2)
	else:
		self_overlap_dict[r.S1]=[r.S2]
df['weighted_distance'] = df.N_mis + df.N_ins*3 + df.N_del*3
df = df.sort_values(['weighted_distance',"total_guideseq_reads"],ascending=[True,False])
def self_overlap_check(s,seq_list,myDict):
	"""
	myDict is a seq:[seq1,seq2,...]

	that means all these seq are self overlapping

	"""

	for i in seq_list:
		if i in myDict:
			if s in myDict[i]:
				return True
	return False
selected_site=[]
for site in df.index.tolist():
	if not self_overlap_check(site,selected_site,self_overlap_dict):
		selected_site.append(site)
# merge the count, optimized
numeric_cols=[ 'total_guideseq_reads', 'total_control_reads', 'total_guideseq_non_dsODN_reads', 'total_control_non_dsODN_reads']+[ 'nuclease.rev.plus', 'nuclease.fwd.plus', 'nuclease.fwd.minus', 'nuclease.rev.minus', 'control.rev.plus', 'control.fwd.plus', 'control.fwd.minus', 'control.rev.minus', 'nuclease.mispriming.plus', 'nuclease.mispriming.minus', 'control.mispriming.plus', 'control.mispriming.minus']
merged_rows = []

for site in selected_site:
	current_site = df.loc[site].copy()
	try:
		others=self_overlap_dict[site]
	except:
		merged_rows.append(current_site)
		continue
	group = df.loc[others+[site]]
	summed = group[numeric_cols].sum()
	current_site[numeric_cols] = summed
	joined = ", ".join(group['merged_start_pos'].astype(str))
	current_site['merged_start_pos'] = joined
	merged_rows.append(current_site)
	# for c in :
		# for i in others:
			# df.at[site,c]+=df.at[i,c]
	# c="merged_start_pos"
	# for i in others:
		# df.at[site,c]+=", "+df.at[i,c]
# df = df.loc[selected_site]
df = pd.DataFrame(merged_rows)
df.to_csv(f"{label}.identified.final.tsv",sep="\t",index=False)

os.system(f"rm {tmp_file}")
os.system(f"rm {tmp_out}")