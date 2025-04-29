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
warnings.filterwarnings("ignore")

def merge_sites(d):
	use_cols=["Site_SubstitutionsOnly.chr","Site_SubstitutionsOnly.start","Site_SubstitutionsOnly.end","Site_SubstitutionsOnly.site_sequence","Site_SubstitutionsOnly.on_target_sequence","Site_SubstitutionsOnly.strand","Site_SubstitutionsOnly.distance","Site_SubstitutionsOnly.cas9_cut_position","Site_SubstitutionsOnly.distance_between_cas9_cut_and_most_frequent_guideseq_read_start","Site_GapsAllowed.chr","Site_GapsAllowed.start","Site_GapsAllowed.end","Site_GapsAllowed.site_sequence","Site_GapsAllowed.on_target_sequence","Site_GapsAllowed.strand","Site_GapsAllowed.distance","Site_GapsAllowed.cas9_cut_position","Site_GapsAllowed.distance_between_cas9_cut_and_most_frequent_guideseq_read_start","total_guideseq_reads","total_control_reads","total_guideseq_non_dsODN_reads","total_control_non_dsODN_reads","total.nuclease.rev.plus","total.nuclease.fwd.plus","total.nuclease.fwd.minus","total.nuclease.rev.minus","total.control.rev.plus","total.control.fwd.plus","total.control.fwd.minus","total.control.rev.minus","total.nuclease.mispriming.plus","total.nuclease.mispriming.minus","total.control.mispriming.plus","total.control.mispriming.minus","merged_guide_seq_read_start_positions_1_based","on_target_sequence","min_cas9_cut_distance","Number_assignments"]
	d=d.sort_values(['min_cas9_cut_distance'])
	dp = d.head(n=1)
	for c in [ 'total_guideseq_reads', 'total_control_reads', 'total_guideseq_non_dsODN_reads', 'total_control_non_dsODN_reads']:
		dp.loc[:, c] = d[c].sum()
	for c in [ 'nuclease.rev.plus', 'nuclease.fwd.plus', 'nuclease.fwd.minus', 'nuclease.rev.minus', 'control.rev.plus', 'control.fwd.plus', 'control.fwd.minus', 'control.rev.minus', 'nuclease.mispriming.plus', 'nuclease.mispriming.minus', 'control.mispriming.plus', 'control.mispriming.minus']:
		dp.loc[:,f"total.{c}"] = d[c].sum()
	c="guide_seq_read_start_positions_1_based"
	dp.loc[:,f'merged_{c}'] = ",".join(d[c].astype(str).tolist())
	dp.loc[:,"Number_assignments"] = d['Number_assignments'].max()
	return dp[use_cols]
label = sys.argv[1]
ncore = int(sys.argv[2])
min_read_count = int(sys.argv[3])
df = pd.read_csv(f"{label}.identified.raw.tsv",sep="\t")
group_sizes = df.groupby(['#chr', 'start', 'end']).size()

# Map the group sizes back to the original DataFrame
df['Number_assignments'] = df.set_index(['#chr', 'start', 'end']).index.map(group_sizes)

df['site_chr'] = df['#chr']
df['site_start'] = np.where(df['Site_SubstitutionsOnly.start'] != "", 
                       df['Site_SubstitutionsOnly.start'], 
                       df['Site_GapsAllowed.start'])
df['site_end'] = np.where(df['Site_SubstitutionsOnly.end'] != "", 
                       df['Site_SubstitutionsOnly.end'], 
                       df['Site_GapsAllowed.end'])
df['site_strand'] = np.where(df['Site_SubstitutionsOnly.strand'] != "", 
                       df['Site_SubstitutionsOnly.strand'], 
                       df['Site_GapsAllowed.strand'])
import time
# start_time = time.time()

df_list = [group for _, group in df.groupby(['site_chr','site_start','site_end','site_strand','on_target_sequence'])]
# elapsed_time = time.time() - start_time

# print ("done groupby") 2.5 min for 9 million sites
# hours = int(elapsed_time // 3600)
# minutes = int((elapsed_time % 3600) // 60)
# seconds = int(elapsed_time % 60)
# print(f"Running time: {hours}h {minutes}m {seconds}s")

out_list = Parallel(n_jobs=ncore,verbose=10,batch_size=10)(delayed(merge_sites)(d) for d in df_list)
df = pd.concat(out_list)

for c in ["Site_SubstitutionsOnly.start","Site_SubstitutionsOnly.end","Site_SubstitutionsOnly.distance","Site_SubstitutionsOnly.cas9_cut_position","Site_SubstitutionsOnly.distance_between_cas9_cut_and_most_frequent_guideseq_read_start","Site_GapsAllowed.start","Site_GapsAllowed.end","Site_GapsAllowed.distance","Site_GapsAllowed.cas9_cut_position","Site_GapsAllowed.distance_between_cas9_cut_and_most_frequent_guideseq_read_start"]:
	df[c] = df[c].replace("",np.nan)
	df[c] = df[c].astype("Int64")

if min_read_count>0:
	df[df.total_guideseq_reads<min_read_count].to_csv(f"{label}.identified.min_read_count.filtered.tsv",sep="\t",index=False)
	df[df.total_guideseq_reads>=min_read_count].to_csv(f"{label}.identified.final.tsv",sep="\t",index=False)
else:
	df.to_csv(f"{label}.identified.final.tsv",sep="\t",index=False)

