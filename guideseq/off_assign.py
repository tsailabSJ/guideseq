import pandas as pd
import numpy as np
import glob
from joblib import Parallel, delayed
import re
import regex
import sys
import gc

def off_target_assignment(df,target_list,flank):
	output_list = []
	for i,r in df.iterrows():
		# left extend sequence
		match = re.match(r"(chr\w+):(\d+)-(\d+)", r["left_extend_coord"])
		chromosome, start_left, end = match.groups()
		start_left = int(start_left)
		# right extend sequence
		match = re.match(r"(chr\w+):(\d+)-(\d+)", r["right_extend_coord"])
		chromosome, start_right, end = match.groups()
		start_right = int(start_right)    
		site_list = r.tolist()
		distance_list = []
		for t_seq in target_list:
			flag = True
			for i in [flank,0]:
				start_i_left = start_left+i
				w_seq_left = r["left_extend_seq"]
				if i!=0:
					w_seq_left = r["left_extend_seq"][i:-i]
				D_left,align1_left,align2_left = target_window_match(w_seq_left, start_i_left,t_seq,max_score=7)
				start_i_right = start_right+i
				w_seq_right = r["right_extend_seq"]
				if i!=0:
					w_seq_right = r["right_extend_seq"][i:-i]
				D_right,align1_right,align2_right = target_window_match(w_seq_right, start_i_right,t_seq,max_score=7)
				if D_left!=100:
					flag = False
					distance_list.append(D_left)
					site_list.append([D_left,align1_left,align2_left])
					# print (t_seq,chromosome,D_left,align1_left,align2_left)
					break
				if D_right!=100:
					flag = False
					distance_list.append(D_right)
					site_list.append([D_right,align1_right,align2_right])
					# print (t_seq,chromosome,D_right,align1_right,align2_right)
					break                
			if flag:
				distance_list.append(100)
				site_list.append([100,[],[]])

		min_distance = min(distance_list)
		try:
			second_min_D = min([x for x in distance_list if x != min_distance])
		except:
			second_min_D=100
		argmin_indices = [i for i, x in enumerate(distance_list) if x == min_distance]
		target_min_list = [target_list[i] for i in argmin_indices]
		if min_distance == 100:
			target_min_list=[]
		site_list.append(min_distance)
		site_list.append(second_min_D)
		site_list.append(target_min_list)
		output_list.append(site_list)   
		del site_list
		gc.collect()
	out = pd.DataFrame(output_list,columns = df.columns.tolist()+target_list+['min_distance','second_min_D','min_distance_target_seq_list'])        
	del output_list
	gc.collect()
	return out

def reverseComplement(sequence):
	if sys.version_info[0] < 3:
		tab = string.maketrans("ACGTacgt", "TGCATGCA")
	else:
		tab = str.maketrans("ACGTacgt", "TGCATGCA")
	return sequence.translate(tab)[::-1]
import regex


def target_window_match(window_sequence, window_start,target_sequence,max_score=7):
			   
	offtarget_sequence_no_bulge, number_mismatch, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
	bulged_offtarget_sequence, length, distance, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end, \
	realigned_target_sequence = alignSequences(target_sequence, window_sequence, max_score)
	
	if chosen_alignment_strand_m == "+":
		non_bulged_target_start_absolute = start_no_bulge + window_start
		non_bulged_target_end_absolute = end_no_bulge + window_start
		non_bulged_cas9_cut_position = non_bulged_target_end_absolute-5
		distance_m = number_mismatch
	elif chosen_alignment_strand_m == "-":
		non_bulged_target_start_absolute = window_start - end_no_bulge+len(window_sequence)
		non_bulged_target_end_absolute = window_start - start_no_bulge+len(window_sequence)
		non_bulged_cas9_cut_position = non_bulged_target_start_absolute+6
		distance_m = number_mismatch
	else:
		non_bulged_target_start_absolute, non_bulged_target_end_absolute,non_bulged_cas9_cut_position = '', '', ''
		distance_m = 100
	if chosen_alignment_strand_b == "+":
		bulged_target_start_absolute = bulged_start + window_start
		bulged_target_end_absolute = bulged_end + window_start
		bulged_cas9_cut_position = bulged_target_end_absolute-5
		distance_b = substitutions+insertions+deletions	
	elif chosen_alignment_strand_b == "-":
		bulged_target_start_absolute = window_start - bulged_end+len(window_sequence)
		bulged_target_end_absolute = window_start - bulged_start+len(window_sequence)
		bulged_cas9_cut_position = bulged_target_start_absolute+6
		distance_b = substitutions+insertions+deletions	
	else:
		bulged_target_start_absolute, bulged_target_end_absolute,bulged_cas9_cut_position = '', '', ''
		distance_b = 100	

	# non bulged alignments
	non_bulged_alignment=[non_bulged_target_start_absolute,non_bulged_target_end_absolute,number_mismatch,offtarget_sequence_no_bulge,chosen_alignment_strand_m,non_bulged_cas9_cut_position]
	
	# bulged alignments	
	bulged_alignment=[bulged_target_start_absolute,bulged_target_end_absolute,substitutions+insertions+deletions,bulged_offtarget_sequence,chosen_alignment_strand_b,substitutions, insertions, deletions,realigned_target_sequence,bulged_cas9_cut_position]	
	
	return min(distance_b,distance_m),non_bulged_alignment,bulged_alignment


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

file = sys.argv[1]
ncore=int(sys.argv[2])
target=sys.argv[3]
target_list = pd.read_csv(target,header=None)[0].tolist()
df = pd.read_csv(file,sep="\t",header=None)
df.columns = ['#chr','start','end',"nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus","control.rev.plus","control.fwd.plus","control.fwd.minus","control.rev.minus","nuclease.mispriming.plus","nuclease.mispriming.minus","control.mispriming.plus","control.mispriming.minus",'left_extend_coord','left_extend_seq','right_extend_coord','right_extend_seq']
# print (df)
# print(df.loc[1])

dfs = np.array_split(df, ncore)
del df
gc.collect()
flank = 10
df_list = Parallel(n_jobs=ncore,verbose=10)(delayed(off_target_assignment)(dfs[i],target_list,flank) for i in range(ncore))
df = pd.concat(df_list)
df.to_csv(f"{file}.off_assigned.tsv",sep="\t",index=False)

