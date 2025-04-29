# IdentifyOffTargetSiteSequences.py
#
# 11-3-2024, rewrite identify code

from __future__ import print_function
import subprocess
import datetime

import argparse
import collections
import numpy
import os
import string
import operator
import pyfaidx
import re
import regex
import logging
from Levenshtein import distance
import pandas as pd
# import dill
logger = logging.getLogger('root')
# from skbio.alignment import StripedSmithWaterman
import numpy as np
# from scipy.signal import argrelextrema
from joblib import Parallel, delayed
logger.propagate = False
import glob
import sys
import os
import time

#-------- HemTools job submission code-----------------------
def submit_job(command):
	print (command)
	pipe = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
	job_id_regex = b"Job <([0-9]+)>"
	output = pipe.communicate()[0]
	match = re.match(job_id_regex,output)
	submission_id = match.group(1).decode('utf-8')
	return submission_id
def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()
def multireplace(myString, myDict):
	## keywords format is {{string}}
	for k in myDict:
		if myDict[k] == None:
			value = ""
		else:
			value = myDict[k]
		myString = str(myString).replace("{{"+str(k)+"}}",str(value))
	return myString


def regexFromSequence(seq, lookahead=True, indels=1, errors=7):
	seq = seq.upper()
	"""
	Given a sequence with ambiguous base characters, returns a regex that matches for
	the explicit (unambiguous) base characters
	"""
	from Bio.Data import IUPACData
	# IUPAC_notation_regex = {'N': '[ATCGN]',
							# 'Y': '[CTY]',
							# 'R': '[AGR]',
							# 'W': '[ATW]',
							# 'S': '[CGS]',
							# 'A': 'A',
							# 'T': 'T',
							# 'C': 'C',
							# 'G': 'G'}
	IUPAC_notation_regex = IUPACData.ambiguous_dna_values
	pattern = ''

	for c in seq:
		pattern += "["+IUPAC_notation_regex[c]+"]"

	if lookahead:
		pattern = '(?b:' + pattern + ')'

	pattern_standard = pattern + '{{s<={0}}}'.format(errors)
	pattern_gap = pattern + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)
	return pattern_standard, pattern_gap

"""
Realigned TargetSequence and OffTargetSequence when indels are present in the local alignment
"""
def extendedPattern(seq, indels=1, errors=7):
	IUPAC_notation_regex_extended = {'N': '[ATCGN]','-': '[ATCGN]','Y': '[CTY]','R': '[AGR]','W': '[ATW]','S': '[CGS]','A': 'A','T': 'T','C': 'C','G': 'G'}
	realign_pattern = ''
	for c in seq:
		realign_pattern += IUPAC_notation_regex_extended[c]
	return '(?b:' + realign_pattern + ')' + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)


def realignedSequences(targetsite_sequence, chosen_alignment, errors=7):
	match_sequence = chosen_alignment.group()
	substitutions, insertions, deletions = chosen_alignment.fuzzy_counts

	# get the .fuzzy_counts associated to the matching sequence after adjusting for indels, where 0 <= INS, DEL <= 1
	realigned_fuzzy = (substitutions, max(0, insertions - 1), max(0, deletions - 1))

	if insertions:  # DNA-bulge
		if targetsite_sequence.index('N') > len(targetsite_sequence)/2:  # PAM is on the right end
			targetsite_realignments = [targetsite_sequence[:i + 1] + '-' + targetsite_sequence[i + 1:] for i in range(targetsite_sequence.index('N') + 1)]
		else:
			targetsite_realignments = [targetsite_sequence[:i] + '-' + targetsite_sequence[i:] for i in range(targetsite_sequence.index('N'), len(targetsite_sequence))]
	else:
		targetsite_realignments = [targetsite_sequence]

	realigned_target_sequence, realigned_offtarget_sequence = None, ''  # in case the matching sequence is not founded

	for seq in targetsite_realignments:
		if deletions:  # RNA-bulge
			match_realignments = [match_sequence[:i + 1] + '-' + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
			match_pattern = [match_sequence[:i + 1] + seq[i + 1] + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
		else:
			match_realignments = match_pattern = [match_sequence]

		x = extendedPattern(seq, errors)
		for y_pattern, y_alignment in zip(match_pattern, match_realignments):
			m = regex.search(x, y_pattern, regex.BESTMATCH)
			if m and m.fuzzy_counts == realigned_fuzzy:
				realigned_target_sequence, realigned_offtarget_sequence = seq, y_alignment
	return realigned_target_sequence, realigned_offtarget_sequence


"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match.
"""
def alignSequences(targetsite_sequence, window_sequence, max_score=7):

	window_sequence = window_sequence.upper()
	query_regex_standard, query_regex_gap = regexFromSequence(targetsite_sequence, errors=max_score)

	# Try both strands
	alignments_mm, alignments_bulge = list(), list()
	alignments_mm.append(('+', 'standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
	alignments_mm.append(('-', 'standard', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
	alignments_bulge.append(('+', 'gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
	alignments_bulge.append(('-', 'gapped', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))

	lowest_distance_score, lowest_mismatch = 100, max_score + 1
	chosen_alignment_b, chosen_alignment_m, chosen_alignment_strand_b, chosen_alignment_strand_m = None, None, '', ''

	for aln_m in alignments_mm:
		strand_m, alignment_type_m, match_m = aln_m
		# print (aln_m)
		if match_m != None:
			mismatches, insertions, deletions = match_m.fuzzy_counts
			if mismatches < lowest_mismatch:
				chosen_alignment_m = match_m
				chosen_alignment_strand_m = strand_m
				lowest_mismatch = mismatches

	for aln_b in alignments_bulge:
		strand_b, alignment_type_b, match_b = aln_b
		if match_b != None:
			substitutions, insertions, deletions = match_b.fuzzy_counts
			if insertions or deletions:
				distance_score = substitutions + (insertions + deletions) * 3
				edistance = substitutions + insertions + deletions
				if distance_score < lowest_distance_score and edistance < lowest_mismatch:
					chosen_alignment_b = match_b
					chosen_alignment_strand_b = strand_b
					lowest_distance_score = distance_score

	if chosen_alignment_m:
		offtarget_sequence_no_bulge = chosen_alignment_m.group()
		mismatches = chosen_alignment_m.fuzzy_counts[0]
		start_no_bulge = chosen_alignment_m.start()
		end_no_bulge = chosen_alignment_m.end()
	else:
		offtarget_sequence_no_bulge, mismatches, start_no_bulge, end_no_bulge, chosen_alignment_strand_m = '', '', '', '', ''

	bulged_offtarget_sequence, score, length, substitutions, insertions, deletions, bulged_start, bulged_end, realigned_target = \
		'', '', '', '', '', '', '', '', 'none'
	if chosen_alignment_b:
		realigned_target, bulged_offtarget_sequence = realignedSequences(targetsite_sequence, chosen_alignment_b, max_score)
		if bulged_offtarget_sequence:
			length = len(chosen_alignment_b.group())
			substitutions, insertions, deletions = chosen_alignment_b.fuzzy_counts
			score = substitutions + (insertions + deletions) * 3
			bulged_start = chosen_alignment_b.start()
			bulged_end = chosen_alignment_b.end()
		else:
			chosen_alignment_strand_b = ''

	return [offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m, start_no_bulge, end_no_bulge,
			bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end, realigned_target]

def hamming_distance(s1, s2):
	if len(s1) != len(s2):
		raise ValueError("Strand lengths are not equal!")
	return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))
def reverseComplement(sequence):
	if sys.version_info[0] < 3:
		tab = string.maketrans("ACGTacgt", "TGCATGCA")
	else:
		tab = str.maketrans("ACGTacgt", "TGCATGCA")
	return sequence.translate(tab)[::-1]

"""
rewrite 11-3-2024
# R2 bed
bedtools_command = f'{bedtools} bamtobed -i {output_UMIdedup_bam_file} | grep "/2" > $label.dedup.R2.bed'
"""
def identify_sites(ref,nuc_bam,con_bam,nuc_mispriming_bam,con_mispriming_bam,label,outdir,target,njobs=None,extend_spacer=17,extend_pam=6,flank=10,bedtools="bedtools",chrom_size=None,min_read_count=0,max_control_reads=5):


	if not os.path.exists(outdir):
		try:
			os.makedirs(outdir)
		except:
			print ("Folder exist")

	# get R2 bed
	bedtools_command = f'{bedtools} bamtobed -i {nuc_bam} | grep "/2" > {outdir}/{label}.nuclease.R2.bed'
	subprocess.call(bedtools_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	bedtools_command = f'{bedtools} bamtobed -i {con_bam} | grep "/2" > {outdir}/{label}.control.R2.bed'
	subprocess.call(bedtools_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	
	bedtools_command = f'{bedtools} bamtobed -i {nuc_mispriming_bam} | grep "/2" > {outdir}/{label}.nuclease.mispriming.R2.bed'
	subprocess.call(bedtools_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	bedtools_command = f'{bedtools} bamtobed -i {con_mispriming_bam} | grep "/2" > {outdir}/{label}.control.mispriming.R2.bed'
	subprocess.call(bedtools_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	
	# summarize R2 start
	# List of bed files
	bed_files = [
		f"{outdir}/{label}.nuclease.R2.bed",
		f"{outdir}/{label}.control.R2.bed"
	]
	bed_files_mispriming = [
		f"{outdir}/{label}.nuclease.mispriming.R2.bed",
		f"{outdir}/{label}.control.mispriming.R2.bed"
	]
	# Commands template for each bed file
	commands = []
	for bed_file in bed_files:
		commands.append(f"awk -F'\\t' '$4 ~ /_rev_/ && $6==\"+\" {{ print $1\"\\t\"$2\"\\t\"$2+1\"\\t1\"}}' {bed_file} | sort -k1,1 -k2,2n | {bedtools} groupby -g 1,2,3 -c 4 -o sum > {bed_file.replace('.R2.bed', '.R2.rev.plus.bdg')}")
		commands.append(f"awk -F'\\t' '$4 ~ /_fwd_/ && $6==\"+\" {{ print $1\"\\t\"$2\"\\t\"$2+1\"\\t1\"}}' {bed_file} | sort -k1,1 -k2,2n | {bedtools} groupby -g 1,2,3 -c 4 -o sum > {bed_file.replace('.R2.bed', '.R2.fwd.plus.bdg')}")
		commands.append(f"awk -F'\\t' '$4 ~ /_fwd_/ && $6==\"-\" {{ print $1\"\\t\"$3-1\"\\t\"$3\"\\t1\"}}' {bed_file} | sort -k1,1 -k2,2n | {bedtools} groupby -g 1,2,3 -c 4 -o sum > {bed_file.replace('.R2.bed', '.R2.fwd.minus.bdg')}")
		commands.append(f"awk -F'\\t' '$4 ~ /_rev_/ && $6==\"-\" {{ print $1\"\\t\"$3-1\"\\t\"$3\"\\t1\"}}' {bed_file} | sort -k1,1 -k2,2n | {bedtools} groupby -g 1,2,3 -c 4 -o sum > {bed_file.replace('.R2.bed', '.R2.rev.minus.bdg')}")
	for bed_file in bed_files_mispriming:
		commands.append(f"awk -F'\\t' '$4 ~ /_no_adapter_/ && $6==\"-\" {{ print $1\"\\t\"$3-1\"\\t\"$3\"\\t1\"}}' {bed_file} | sort -k1,1 -k2,2n | {bedtools} groupby -g 1,2,3 -c 4 -o sum > {bed_file.replace('.R2.bed', '.R2.minus.bdg')}")
		commands.append(f"awk -F'\\t' '$4 ~ /_no_adapter_/ && $6==\"+\" {{ print $1\"\\t\"$2\"\\t\"$2+1\"\\t1\"}}' {bed_file} | sort -k1,1 -k2,2n | {bedtools} groupby -g 1,2,3 -c 4 -o sum > {bed_file.replace('.R2.bed', '.R2.plus.bdg')}")
	Parallel(n_jobs=njobs,verbose=10)(delayed(os.system)(c) for c in commands)
	
	# merge output
	union_bedg = f"{bedtools} unionbedg -i {outdir}/{label}.nuclease.R2.rev.plus.bdg {outdir}/{label}.nuclease.R2.fwd.plus.bdg {outdir}/{label}.nuclease.R2.fwd.minus.bdg {outdir}/{label}.nuclease.R2.rev.minus.bdg {outdir}/{label}.control.R2.rev.plus.bdg {outdir}/{label}.control.R2.fwd.plus.bdg {outdir}/{label}.control.R2.fwd.minus.bdg {outdir}/{label}.control.R2.rev.minus.bdg {outdir}/{label}.nuclease.mispriming.R2.plus.bdg {outdir}/{label}.nuclease.mispriming.R2.minus.bdg {outdir}/{label}.control.mispriming.R2.plus.bdg {outdir}/{label}.control.mispriming.R2.minus.bdg > {outdir}/{label}.R2.count.tsv"
	logger.info(f'Running union_bedg command: {union_bedg}' )
	subprocess.call(union_bedg, shell=True,stdout=sys.stdout,stderr=sys.stderr)

	# getfasta
	# seq_extend = f"awk '{{print $1, $2 - {extend_spacer} - {flank}, $3 + {extend_pam} + {flank}}}' OFS='\t' {outdir}/{label}.R2.count.tsv > {outdir}/{label}.seq1.bed"
	seq_extend = f"awk '{{start = ($2 - {extend_spacer} - {flank} < 0) ? 0 : $2 - {extend_spacer} - {flank}; print $1, start, $3 + {extend_pam} + {flank}}}' OFS='\t' {outdir}/{label}.R2.count.tsv | {bedtools} intersect -a - -b {chrom_size} > {outdir}/{label}.seq1.bed"
	logger.info(f'Running seq_extend command: {seq_extend}' )
	subprocess.call(seq_extend, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	getseq = f"{bedtools} getfasta -fi {ref} -bed {outdir}/{label}.seq1.bed -fo {outdir}/{label}.seq1.tsv -tab"
	subprocess.call(getseq, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	# seq_extend = f"awk '{{print $1, $2 - {extend_pam} - {flank}, $3 + {extend_spacer} + {flank}}}' OFS='\t' {outdir}/{label}.R2.count.tsv > {outdir}/{label}.seq2.bed"
	seq_extend = f"awk '{{start = ($2 - {extend_pam} - {flank} < 0) ? 0 : $2 - {extend_pam} - {flank}; print $1, start, $3 + {extend_spacer} + {flank}}}' OFS='\t' {outdir}/{label}.R2.count.tsv | {bedtools} intersect -a - -b {chrom_size}  > {outdir}/{label}.seq2.bed"
	logger.info(f'Running seq_extend command: {seq_extend}' )
	subprocess.call(seq_extend, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	getseq = f"{bedtools} getfasta -fi {ref} -bed {outdir}/{label}.seq2.bed -fo {outdir}/{label}.seq2.tsv -tab"
	subprocess.call(getseq, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	
	# append to merged output
	command = f"paste {outdir}/{label}.R2.count.tsv {outdir}/{label}.seq1.tsv {outdir}/{label}.seq2.tsv > {outdir}/{label}.raw.tsv"
	subprocess.call(command, shell=True,stdout=sys.stdout,stderr=sys.stderr)


	# split raw.tsv, run about 1 day per 1M data, 10 cores, low mem, high-cpu job
	logger.info(f'Running speed_up mode' )
	N=200000
	command = f"mkdir -p {outdir}/speed_up_{label};split -l {N} {outdir}/{label}.raw.tsv {outdir}/speed_up_{label}/{label}.raw.split_"
	# print (command)
	subprocess.call(command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	ncore=10
	current_script_path = os.path.abspath(__file__)
	python = "/research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/anaconda3/envs/changeseq_py3/bin/python"
	myDate=str(datetime.date.today())
	os.system(f"mkdir -p identify_parallel_log_{myDate}")
	# command = 'bsub -R "rusage[mem=3000] span[hosts=1]" -n 10 -P {{label}} -J {{label}} -q standard -oo identify_parallel_log_{{myDate}}/{{label}}.log -eo identify_parallel_log_{{myDate}}/{{label}}.err {{python}} {{myscript}} {{args}}'
	command = 'bsub -R "rusage[mem=1000] span[hosts=1]" -n {{ncore}} -P {{label}} -J {{label}} -q standard -oo identify_parallel_log_{{myDate}}/{{label}}.log -eo identify_parallel_log_{{myDate}}/{{label}}.{{label2}}.err {{python}} {{myscript}} {{args}}'
	myDict={}
	myDict['python'] = python
	myDict['ncore'] = ncore
	myDict['myDate'] = myDate
	# off_identify.py
	myDict['label']="off_identify_v3"
	myscript = current_script_path.replace("identifyOfftargetSites.py",f"{myDict['label']}.py")
	myDict['myscript']=myscript
	jobIDs = []
	files = glob.glob(f"{outdir}/speed_up_{label}/{label}.raw.split_*")
	for f in files:
		myDict['args']=f"{f} {ncore} {target}"
		myDict['label2']=f.split("_")[-1]
		jobIDs.append(submit_job(multireplace(command,myDict)))
	monitor_jobs(jobIDs)
	
	# merge idenfied
	merge_command = "{ head -n 1 $(ls %s/speed_up_%s/*.identified.tsv | head -n 1); tail -n +2 -q %s/speed_up_%s/*.identified.tsv; } > %s/%s.identified.raw.tsv"%(outdir,label,outdir,label,outdir,label)
	os.system(merge_command)
	# merge unidentified
	merge_command = "{ head -n 1 $(ls %s/speed_up_%s/*.unidentified.tsv | head -n 1); tail -n +2 -q %s/speed_up_%s/*.unidentified.tsv; } > %s/%s.unidentified.tsv"%(outdir,label,outdir,label,outdir,label)
	os.system(merge_command)
	# merge noGUIDEseqReads
	merge_command = "{ head -n 1 $(ls %s/speed_up_%s/*.noGUIDEseqReads.tsv | head -n 1); tail -n +2 -q %s/speed_up_%s/*.noGUIDEseqReads.tsv; } > %s/%s.noGUIDEseqReads.tsv"%(outdir,label,outdir,label,outdir,label)
	os.system(merge_command)

	# off_count_merge.py # fast high-mem job
	ncore=2
	current_script_path = os.path.abspath(__file__)
	python = "/research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/anaconda3/envs/changeseq_py3/bin/python"
	myDate=str(datetime.date.today())
	os.system(f"mkdir -p identify_parallel_log_{myDate}")
	command = 'bsub -R "rusage[mem={{mem}}] span[hosts=1]" -n {{ncore}} -P {{label}} -J {{label}} -q standard -oo identify_parallel_log_{{myDate}}/{{label}}.log -eo identify_parallel_log_{{myDate}}/{{label}}.{{label2}}.err {{python}} {{myscript}} {{args}}'
	myDict={}
	myDict['python'] = python
	myDict['myDate'] = myDate
	myDict['ncore'] = ncore
	
	df = pd.read_csv(f"{outdir}/{label}.identified.raw.tsv",sep="\t")
	myDict['label']="off_count_merge_v2"
	myscript = current_script_path.replace("identifyOfftargetSites.py",f"{myDict['label']}.py")
	myDict['myscript']=myscript
	jobIDs = []
	for s,d in df.groupby("#chr"):
		raw_tsv = f"{outdir}/{label}.{s}.identified.raw.tsv"
		d.to_csv(raw_tsv,sep="\t",index=False)
		myDict['args']=f"{outdir}/{label}.{s}.identified.raw.tsv {ncore}"
		myDict['label2']=s
		if s in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX","chr8","chr9","chr11","chr10"]:
			myDict['mem'] = 60000
		else:
			myDict['mem'] = 30000
		jobIDs.append(submit_job(multireplace(command,myDict)))
	monitor_jobs(jobIDs)


	get_final_identified(outdir,label)
	for s,d in df.groupby("#chr"):
		raw_tsv = f"{outdir}/{label}.{s}.identified.raw.tsv"
		os.system(f"rm {raw_tsv}")
		raw_tsv = f"{outdir}/{label}.{s}.identified.final.tsv"
		os.system(f"rm {raw_tsv}")
	
	
	
def get_final_identified(outdir,label):

	files = glob.glob(f"{outdir}/{label}.*.identified.final.tsv")
	df = pd.concat([pd.read_csv(f,sep="\t") for f in files])
	Gini_cols=['total.nuclease.rev.plus','total.nuclease.fwd.plus','total.nuclease.fwd.minus','total.nuclease.rev.minus']
	df['dsODN_read_distribution_Gini_index'] = df[Gini_cols].swifter.progress_bar(False).apply(gini_index,axis=1).fillna(1)
	df.start = df.start_absolute
	df.end = df.end_absolute
	df = df.drop(['start_absolute','end_absolute'],axis=1)
	df['site_name'] = df['#chr']+":"+df.start.astype(str)+"-"+df.end.astype(str)+"_"+df.on_target_sequence
	df.to_csv(f"{outdir}/{label}.identified.unfiltered.final.tsv",sep="\t",index=False)
	df_hc = get_hc_sites(df)
	df_hc.to_csv(f"{outdir}/{label}.matched.final.high_confidence.tsv",index=False,sep="\t")

def get_hc_sites(df):
	Gini_cols=['total.nuclease.rev.plus','total.nuclease.fwd.plus','total.nuclease.fwd.minus','total.nuclease.rev.minus']
	# bi-direction evidence
	df = df[df[Gini_cols].astype(bool).astype(int).sum(axis=1)>1]
	site_names1 = df[df.editDistance<=4].site_name.tolist()
	site_names2 = df[(df.editDistance>4)&(df.total_guideseq_reads>=10)&(df.distance_between_cas9_cut_and_most_frequent_guideseq_read_start<=10)&(df.dsODN_read_distribution_Gini_index<0.5)].site_name.tolist()
	return df[df.site_name.isin(site_names1+site_names2)]


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

	
def monitor_jobs(job_ids):
	"""Monitor a list of job IDs until all jobs are finished."""
	while job_ids:
		for jobid in job_ids[:]:  # Iterate over a copy of the list
			if is_job_finished(jobid):
				print(f"Job {jobid} has finished.")
				job_ids.remove(jobid)  # Remove finished job from the list
		if job_ids:
			print(f"Waiting for {len(job_ids)} job(s) to finish: {job_ids}")
			time.sleep(100)  # Wait 10 seconds before checking again
	print("All jobs are finished!")

def is_job_finished(jobid):
	try:
		# Run the bjobs command to check job status
		result = subprocess.Popen("bjobs %s"%(str(jobid)), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).communicate()[0]
		result = result.decode('utf-8')
		# If the job is not found, it's likely finished
		if "not found" in result or "EXIT" in result or "DONE" in result:
			return True
		return False
	except Exception as e:
		print(f"Error checking job status: {e}")
		return False


def save_object(obj,out):

	try:
		with open(out, "wb") as f:
			dill.dump(obj, f)
	except Exception as ex:
		print("Error during pickling object (Possibly unsupported):", ex)
 
def py2min(myList):
	out = [i  for i in myList if i != "" ]
	return min(out)
	
def py2max(myList):
	out = [i  for i in myList if i != "" ]
	return max(out)
def revcomp(seq):
	try: ## python2
		tab = string.maketrans(b"ACTG", b"TGAC")
	except:  ## python3
		tab = bytes.maketrans(b"ACTG", b"TGAC")
	return seq.translate(tab)[::-1]

def loadFileIntoArray(filename):
	with open(filename, 'rU') as f:
		keys = f.readline().rstrip('\r\n').split('\t')[1:]
		data = collections.defaultdict(dict)
		for line in f:
			filename, rest = processLine(line)
			line_to_dict = dict(zip(keys, rest))
			data[filename] = line_to_dict
	return data



def main():
	parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
	parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
	parser.add_argument('--nuclease_sam', help='SAM file', required=True)
	parser.add_argument('--control_sam', help='SAM file', required=True)
	parser.add_argument('--label', help='File to output identified sites to.', required=True)
	parser.add_argument('--max_score', help='Score threshold', type=int, default=7)
	parser.add_argument('--target',  help='can be a file or a single str', required=True)

	args = parser.parse_args()

	# identify_sites(args.samfile[0], args.ref, args.outfile, annotations, args.window, args.max_score,args.control_primer)

if __name__ == "__main__":
	main()
