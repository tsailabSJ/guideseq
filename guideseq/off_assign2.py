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

file=sys.argv[1]

df = pd.read_csv(file,sep="\t")
df['min_distance_target_seq_list'] = df['min_distance_target_seq_list'].apply(eval)

# # off-target assignment
# 1. use min_distance< second_min_D
# 2. always allow multiple off-targets in one region, choose the closest cas9 as true off-target, put others into "others" column


total = df.shape[0]
df[df.min_distance>=df.second_min_D][df.columns[:19]].to_csv(f"{file}.unidentified.tsv",sep="\t",index=False)
df = df[df.min_distance<df.second_min_D]
keep = df.shape[0]

print (f"Before filter {total}")
print (f"After filter {keep}, {keep/total}")


# # find true off-target


out = []
for i,r in df.iterrows():
    off_target_list = r['min_distance_target_seq_list']
    current_read_start = r[df.columns[:15]].tolist()
    for off in off_target_list:
        off_row = eval(r[off])
        out.append(current_read_start+off_row[1]+off_row[2]+[off])

df2 = pd.DataFrame(out)

df2.columns = df.columns[:15].tolist()+["Site_SubstitutionsOnly.start","Site_SubstitutionsOnly.end","Site_SubstitutionsOnly.distance","Site_SubstitutionsOnly.site_sequence","Site_SubstitutionsOnly.strand","Site_SubstitutionsOnly.cas9_cut_position","Site_GapsAllowed.start","Site_GapsAllowed.end","Site_GapsAllowed.distance","Site_GapsAllowed.site_sequence","Site_GapsAllowed.strand","Site_GapsAllowed.substitutions","Site_GapsAllowed.insertions","Site_GapsAllowed.deletions","Site_GapsAllowed.on_target_sequence","Site_GapsAllowed.cas9_cut_position","on_target_sequence"]



def get_detailed_info(r):
    r['Site_SubstitutionsOnly.chr'] = r['#chr']
    r['Site_GapsAllowed.chr'] = r['#chr']
    r['Site_SubstitutionsOnly.distance_between_cas9_cut_and_most_frequent_guideseq_read_start']=""
    r['Site_GapsAllowed.distance_between_cas9_cut_and_most_frequent_guideseq_read_start']=""
    if r['Site_SubstitutionsOnly.site_sequence'] !="":
        r['Site_SubstitutionsOnly.distance_between_cas9_cut_and_most_frequent_guideseq_read_start']=abs(r['end']-r['Site_SubstitutionsOnly.cas9_cut_position'])
    if r['Site_GapsAllowed.site_sequence'] !="":
        r['Site_GapsAllowed.distance_between_cas9_cut_and_most_frequent_guideseq_read_start']=abs(r['end']-r['Site_GapsAllowed.cas9_cut_position'])
    r['total_guideseq_reads'] = r[["nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus"]].sum()
    r['total_control_reads'] = r[["control.rev.plus","control.fwd.plus","control.fwd.minus","control.rev.minus"]].sum()
    r['total_guideseq_non_dsODN_reads'] = r[["nuclease.mispriming.plus","nuclease.mispriming.minus"]].sum()
    r['total_control_non_dsODN_reads'] = r[["control.mispriming.plus","control.mispriming.minus"]].sum()
    r['guide_seq_read_start_positions_1_based'] = r['end']
    return r
df3 = df2.swifter.progress_bar(False).apply(get_detailed_info,axis=1)


use_cols=["#chr","start","end","Site_SubstitutionsOnly.chr","Site_SubstitutionsOnly.start","Site_SubstitutionsOnly.end","Site_SubstitutionsOnly.site_sequence","Site_SubstitutionsOnly.on_target_sequence","Site_SubstitutionsOnly.strand","Site_SubstitutionsOnly.distance","Site_SubstitutionsOnly.cas9_cut_position","Site_SubstitutionsOnly.distance_between_cas9_cut_and_most_frequent_guideseq_read_start","Site_GapsAllowed.chr","Site_GapsAllowed.start","Site_GapsAllowed.end","Site_GapsAllowed.site_sequence","Site_GapsAllowed.on_target_sequence","Site_GapsAllowed.strand","Site_GapsAllowed.distance","Site_GapsAllowed.cas9_cut_position","Site_GapsAllowed.distance_between_cas9_cut_and_most_frequent_guideseq_read_start","total_guideseq_reads","total_control_reads","total_guideseq_non_dsODN_reads","total_control_non_dsODN_reads","nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus","control.rev.plus","control.fwd.plus","control.fwd.minus","control.rev.minus","nuclease.mispriming.plus","nuclease.mispriming.minus","control.mispriming.plus","control.mispriming.minus","guide_seq_read_start_positions_1_based","on_target_sequence"]

df3['Site_SubstitutionsOnly.on_target_sequence'] = df3['on_target_sequence']
df3 = df3[use_cols]



df3['min_cas9_cut_distance'] = df3[['Site_SubstitutionsOnly.distance_between_cas9_cut_and_most_frequent_guideseq_read_start','Site_GapsAllowed.distance_between_cas9_cut_and_most_frequent_guideseq_read_start']].replace('',100).min(axis=1)



df4 = df3.sort_values(['#chr','start','min_cas9_cut_distance']).reset_index(drop=True)



df5 = df4.drop_duplicates(['#chr','start'])





df4['off_target_assignment'] = 0
df4.loc[df4.index.isin(df5.index), 'off_target_assignment'] = 1




df4.to_csv(f"{file}.identified.tsv",sep="\t",index=False)
















