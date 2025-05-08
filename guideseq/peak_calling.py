
# python /research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/src/guideseq_pool2/guideseq/peak_calling.py C552.unidentified.tsv cas9
# python /research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/src/guideseq_pool2/guideseq/peak_calling.py C552.noGUIDEseqReads.tsv Control
# python /research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/src/guideseq_pool2/guideseq/peak_calling.py C552.unidentified.tsv cas9 ~/Data/Blacklist/lists/hg38-blacklist.v2.bed.gz ~/Data/Human/hg38/fasta/hg38.fa NGG hg38


import sys
import pandas as pd
import os
import numpy as np
SEACR="/home/yli11/.conda/envs/captureC/bin/SEACR_1.3.sh"


input=sys.argv[1]
myType=sys.argv[2]
BL=sys.argv[3]
fasta=sys.argv[4]
PAM=sys.argv[5]
genome=sys.argv[6]

label = input.split(".")[0]
df = pd.read_csv(input,sep="\t")
G_cols=["nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus"]
C_cols=["control.rev.plus","control.fwd.plus","control.fwd.minus","control.rev.minus"]
df['G_sum'] = df[G_cols].sum(axis=1)
df['C_sum'] = df[C_cols].sum(axis=1)

if myType == "cas9":
	bdg_file=f"{label}.{myType}.bdg"
	peak_file=f"Peak.{label}.{myType}"
	df[['#chr','start','end','G_sum']].sort_values(['#chr','start']).to_csv(bdg_file,sep="\t",header=False,index=False)
	command = f"{SEACR} {bdg_file} 0.01 non stringent {peak_file}"
	os.system(command)
else:
	bdg_file=f"{label}.{myType}.bdg"
	peak_file=f"Peak.{label}.{myType}"
	df = df[df.C_sum>0]
	df[['#chr','start','end','C_sum']].sort_values(['#chr','start']).to_csv(bdg_file,sep="\t",header=False,index=False)
	command = f"{SEACR} {bdg_file} 0.01 non stringent {peak_file}"
	os.system(command)
peak_file += ".stringent.bed"
min_read=10
df2 = pd.read_csv(peak_file,sep="\t",header=None)
df2 = df2[df2[3]>=min_read]
df2.columns = ['#chr','start','end','Total_read','Max_read','Window']
df2.to_csv(peak_file,index=False,sep="\t")
peak_out=f"Peak.{label}.{myType}.rmblck.stringent.bed"
command = f"bedtools intersect -a {peak_file} -b {BL} -v > {peak_out};rm {peak_file}"
os.system(command)
# identify sequence
if myType != "cas9":
	exit()
peak_file = peak_out
peak_out=f"Peak.{label}.{myType}.rmblck.stringent.tsv"
command = f"bedtools intersect -a {peak_file} -b {input} -wo > {peak_out};rm {peak_file}"
os.system(command)
df2 = pd.read_csv(peak_out,sep="\t",header=None,usecols=list(range(26)))
df = df.drop(['G_sum','C_sum'],axis=1)
df2.columns = ['#chr1','start1','end1','Total_read','Max_read','Window']+df.columns.tolist()+['bp']
def gini_index(data):
	# Ensure data is sorted in ascending order
	data = np.sort(np.array(data))
	if sum(data)<=1:
		return 1
	# Calculate the mean of the data
	mean = np.mean(data)
	
	# Total number of observations
	n = len(data)
	
	# Gini Index formula
	gini_sum = sum((i + 1) * data[i] for i in range(n))
	gini = (2 * gini_sum) / (n * sum(data)) - (n + 1) / n
	
	return gini
df2 = df2.groupby(['#chr1','start1','end1','Total_read','Max_read']).sum().reset_index()
df2['dsODN_read_distribution_Gini_index'] = df2[G_cols].apply(gini_index,axis=1).fillna(1)
df2 = df2[df2['dsODN_read_distribution_Gini_index']<0.75]


import pyfaidx
import regex
fasta = pyfaidx.Fasta(fasta)
def reverse_complement(seq):
	complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N', 'R':'Y', 'H':'D'}
	return ''.join(complement[base] for base in reversed(seq))

def regex_from_sequence(seq, insertions=0, deletions=0, mismatches=0):
	pattern = ''.join([{'A':'A', 'T':'T', 'C':'C', 'G':'G', 'N':'[ATCGN]', 'R':'[AG]', 'H':'[ACT]', 'Y':'[TC]', 'D':'[AGT]'}[c] for c in seq.upper()])
	pattern = f'(?b:{pattern})' + f'{{s<={mismatches},i<={insertions},d<={deletions}}}'
	return pattern

def get_genomic_sequence(chromosome, start, end):
	seq = fasta[chromosome][start:end].seq
	return seq

def align_sequences(target_seq, window_seq, window_start,max_insertions=0, max_deletions=0, max_mismatches=0):
	pattern = regex_from_sequence(target_seq, max_insertions, max_deletions, max_mismatches)
	matches = list(regex.finditer(pattern, window_seq))
	return_list = []
	for match in matches:
		start, end = match.span()
		strand = "+"
		num_mismatches, num_insertions, num_deletions = match.fuzzy_counts
		return_list.append([window_seq[start:end], start+window_start, end+window_start, strand])
	
	target_seq_rc = reverse_complement(target_seq)
	pattern = regex_from_sequence(target_seq_rc, max_insertions, max_deletions, max_mismatches)
	matches = list(regex.finditer(pattern, window_seq))
	for match in matches:
		start, end = match.span()
		strand = "-"
		num_mismatches, num_insertions, num_deletions = match.fuzzy_counts
		return_list.append([reverse_complement(window_seq[start:end]), start+window_start, end+window_start, strand, end-start, num_mismatches, num_insertions, num_deletions])
	return return_list #best_match, start, end, strand, length, num_mismatches, num_insertions, num_deletions


def get_off_target(r,PAM):
	chromosome = r['#chr1']
	position = int((r['start1']+r['end1'])/2)
	target_seq = "N"*20+PAM
	window_size=25
	window_start = position - window_size
	window_end = position + window_size
	genomic_sequence = get_genomic_sequence(chromosome, window_start, window_end)
	genomic_sequence = genomic_sequence.upper()
	return_list = align_sequences(target_seq, genomic_sequence, window_start,max_insertions = 0, max_deletions = 0, max_mismatches = 0)
	r['off_target_list'] = return_list
	r['genomic_sequence'] = genomic_sequence
	return r

peak_out=f"Peak.{label}.{myType}.rmblck.stringent.bed"

df2 = df2.apply(lambda r:get_off_target(r,PAM),axis=1)
df2 = df2[['#chr1','start1','end1','off_target_list','Total_read']]
# df2['ID'] = df['#chr1']+":"+df.start1.astype(str)+"-"+df.end1.astype(str)
df2['strand']="+"
df2.columns = ['#chr','start','end','off_target_list','Total_read','strand']
df2.to_csv(peak_out,sep="\t",index=False)

homer_out=f"Peak.{label}.{myType}.rmblck.stringent.annot.tsv"
command = f"annotatePeaks.pl {peak_out} {genome} > {homer_out}"
os.system(command)
# df3 = pd.read_csv(homer_out,sep="\t",index_col=0)
# df3 = df3[df3]
