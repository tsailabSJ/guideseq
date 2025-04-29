import pandas as pd
import numpy as np
import glob
from joblib import Parallel, delayed
import sys
import gc
import swifter
import edlib
import os
import re
import regex
from Bio.Data import IUPACData
def reverseComplement(sequence):
	if sys.version_info[0] < 3:
		tab = string.maketrans("ACGTacgt", "TGCATGCA")
	else:
		tab = str.maketrans("ACGTacgt", "TGCATGCA")
	return sequence.translate(tab)[::-1]
	
# 3/7/2025 updated using edlib
# https://github.com/Martinsos/edlib/issues/88
degenNuc = [("R", "A"), ("R", "G"), 
			("M", "A"), ("M", "C"),
			("W", "A"), ("W", "T"),
			("S", "C"), ("S", "G"),
			("Y", "C"), ("Y", "T"),
			("K", "G"), ("K", "T"),
			("V", "A"), ("V", "C"), ("V", "G"),
			("H", "A"), ("H", "C"), ("H", "T"),
			("D", "A"), ("D", "G"), ("D", "T"),
			("B", "C"), ("B", "G"), ("B", "T"),
			("N", "G"), ("N", "A"), ("N", "T"), ("N", "C"),
			("X", "G"), ("X", "A"), ("X", "T"), ("X", "C")]
def getNiceAlignment2(alignResult, query, target, gapSymbol="-"):
	""" Output alignments from align() in NICE format
	Yichao Li
	add MID counts and aligned start and end position
	
	@param {dictionary} alignResult, output of the method align() 
		NOTE: The method align() requires the argument task="path"
	@param {string} query, the exact query used for alignResult
	@param {string} target, the exact target used for alignResult
	@param {string} gapSymbol, default "-"
		String used to represent gaps in the alignment between query and target
	@return Alignment in NICE format, which is human-readable visual representation of how the query and target align to each other. 
		e.g., for "telephone" and "elephant", it would look like:
		   telephone
			|||||.|.
		   -elephant
		It is represented as dictionary with following fields:
		  - {string} query_aligned
		  - {string} matched_aligned ('|' for match, '.' for mismatch, ' ' for insertion/deletion)
		  - {string} target_aligned
		Normally you will want to print these three in order above joined with newline character.
	"""
	if type(alignResult) is not dict:
		raise Exception("The object alignResult is expected to be a python dictionary. Please check the input alignResult.")

	if 'locations' not in alignResult.keys():
		raise Exception("The object alignResult is expected to contain a field 'locations'. Please check the input alignResult.")

	target_pos = alignResult["locations"][0][0]
	if target_pos == None:
		target_pos = 0
	query_pos = 0 ## 0-indexed
	target_aln = match_aln = query_aln = ""

	if 'cigar' not in alignResult.keys():
		raise Exception("The object alignResult is expected to contain a CIGAR string. Please check the input alignResult.")

	cigar = alignResult["cigar"]

	if alignResult["cigar"] == '' or alignResult["cigar"] == None:
		raise Exception("The object alignResult contains an empty CIGAR string. Users must run align() with task='path'. Please check the input alignResult.")		

	###
	## cigar parsing, motivated by yech1990: https://github.com/Martinsos/edlib/issues/127
	##
	## 'num_occurences' == Number of occurrences of the alignment operation
	## 'alignment_operation' == Cigar symbol/code that represent an alignment operation
	## 
	# 3=1X2=1I1=1I3=2X1=1I1=2I1=
	N_mis=0
	N_ins=0
	N_del=0
	N_match=0
	for num_occurrences, alignment_operation in re.findall("(\d+)(\D)", cigar):
		num_occurrences = int(num_occurrences)
		if alignment_operation == "=":
			N_match+=num_occurrences
			target_aln += target[target_pos : target_pos + num_occurrences]
			target_pos += num_occurrences
			query_aln += query[query_pos : query_pos + num_occurrences]
			query_pos += num_occurrences
			match_aln += "|" * num_occurrences
		elif alignment_operation == "X":
			N_mis+=num_occurrences
			target_aln += target[target_pos : target_pos + num_occurrences]
			target_pos += num_occurrences
			query_aln += query[query_pos : query_pos + num_occurrences]
			query_pos += num_occurrences
			match_aln += "." * num_occurrences
		elif alignment_operation == "D":
			N_del+=num_occurrences
			target_aln += target[target_pos : target_pos + num_occurrences]
			target_pos += num_occurrences
			query_aln += gapSymbol * num_occurrences
			query_pos += 0
			match_aln +=gapSymbol * num_occurrences
		elif alignment_operation == "I":
			N_ins+=num_occurrences
			target_aln += gapSymbol * num_occurrences
			target_pos += 0
			query_aln += query[query_pos : query_pos + num_occurrences]
			query_pos += num_occurrences
			match_aln += gapSymbol * num_occurrences
		else:
			raise Exception("The CIGAR string from alignResult contains a symbol not '=', 'X', 'D', 'I'. Please check the validity of alignResult and alignResult.cigar")

	alignments = {
		'N_mis': N_mis,
		'N_ins': N_ins,
		'N_del': N_del,
		'editDistance': alignResult["editDistance"],
		'start': alignResult["locations"][0][0],
		'end': alignResult["locations"][0][0]+N_match+N_mis+N_del,
		'query_aligned': query_aln,
		'matched_aligned': match_aln,
		'target_aligned': target_aln
	}

	return alignments

def regexFromSequence(seq, lookahead=True,errors=7):
	seq = seq.upper()
	"""
	Given a sequence with ambiguous base characters, returns a regex that matches for
	the explicit (unambiguous) base characters
	"""
	
	IUPAC_notation_regex = IUPACData.ambiguous_dna_values
	pattern = ''

	for c in seq:
		pattern += "["+IUPAC_notation_regex[c]+"]"

	if lookahead:
		pattern = '(?b:' + pattern + ')'

	pattern_standard = pattern + '{{s<={0}}}'.format(errors)
	return pattern_standard
def gapless_alignment(gRNA,window_seq):
	return regex.search(regexFromSequence(gRNA), window_seq, regex.BESTMATCH)


def edlibAlign(gRNA,window_seq,window_seq_revcomp,max_score=7):
	"""align on-target with PAM to window sequence
	
	gRNA: str, 20bp spacer + PAM, e.g., NGG
	
	window_seq, str, this is typically a longer sequence
	
	output
	------
	0 start
	1 end
	2 gRNA_aligned
	3 window_seq_aligned
	4 strand
	5 editDistance
	6 N_mismatch
	7 N_insertion
	8 N_deletion

	
	"""
	output = [None]*9
	result=edlib.align(gRNA,window_seq, mode = "HW", task = "path",k=max_score,additionalEqualities=degenNuc)
	result_revcomp=edlib.align(gRNA,window_seq_revcomp, mode = "HW", task = "path",k=max_score,additionalEqualities=degenNuc)
	if result['editDistance']==-1: # no alignment in + strand
		if result_revcomp['editDistance']==-1: # no alignment given max_score
			return output
		else:
			strand="-"
			alignment=getNiceAlignment2(result_revcomp, gRNA, window_seq_revcomp)
	else: # + strand alignment found
		if result_revcomp['editDistance']==-1: # no alignment in - strand
			strand="+"
			alignment=getNiceAlignment2(result, gRNA, window_seq)
		else: # both aligned, select the best
			if result['editDistance']<=result_revcomp['editDistance']:
				strand="+"
				alignment=getNiceAlignment2(result, gRNA, window_seq)
			else:
				strand="-"
				alignment=getNiceAlignment2(result_revcomp, gRNA, window_seq_revcomp)
	if strand != None:
		N_indel = alignment['N_ins']+alignment['N_del']
		if strand == "+" and N_indel>0:
			# try no gap allowed alignment
			regex_out = gapless_alignment(gRNA,window_seq)
		elif strand == "-" and N_indel>0:
			# try no gap allowed alignment
			regex_out = gapless_alignment(gRNA,window_seq_revcomp)
		else:
			regex_out=None
		if regex_out!=None:
			# gapless alignment found!
			alignment['start'] = regex_out.span()[0]
			alignment['end'] = regex_out.span()[1]
			alignment['query_aligned'] = gRNA
			alignment['target_aligned'] = regex_out.group()
			alignment['editDistance'] = regex_out.fuzzy_counts[0]
			alignment['N_mis'] = regex_out.fuzzy_counts[0]
			alignment['N_ins'] = 0
			alignment['N_del'] = 0
			
		N_indel = alignment['N_ins']+alignment['N_del']
		if N_indel>2: # unlikely tobe true off-target
			strand=None

	output[0] = alignment['start']
	output[1] = alignment['end']
	output[2] = alignment['query_aligned']
	output[3] = alignment['target_aligned']
	output[4] = strand
	output[5] = alignment['editDistance']
	output[6] = alignment['N_mis']
	output[7] = alignment['N_ins']
	output[8] = alignment['N_del']

	
	
	return output

def target_window_match(window_sequence, window_start,target_sequence,max_score=7):
	"""
	window_start, 0-based
	returned start, start_absolute are 0-based, end_absolute is 1-based, all follow bed format
	calculated cas9_cut_position, 1-based, this is the first base on the PAM side after cas9 cut
	"""
	if len(window_sequence)<20:
		return [-1]*10
	wseq_revcomp = reverseComplement(window_sequence)
	
	try:
		start,end,query_aligned,target_aligned,strand,editDistance,N_mis,N_ins,N_del = edlibAlign(target_sequence,window_sequence,wseq_revcomp,max_score=max_score)
	except Exception as e:
		print (e)
		print ("Error",target_sequence, window_sequence, max_score)
		# exit()
	

	if strand == "+":
		start_absolute = start + window_start
		end_absolute = end + window_start
		cas9_cut_position = end_absolute-5
	elif strand == "-":
		start_absolute = window_start - end+len(window_sequence)
		end_absolute = window_start - start+len(window_sequence)
		cas9_cut_position = start_absolute+6
	else:
		editDistance=-1
		start_absolute=-1
		end_absolute=-1
		cas9_cut_position=-1

	out=[start_absolute,end_absolute,query_aligned,target_aligned,strand,editDistance,N_mis,N_ins,N_del,cas9_cut_position]
	
	return out


def off_target_assignment(df,target_list,flank,max_score=7):
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
		off_target_list = [] # identified off-target lst in a list
		distance_list = [] # for each target, we have its distance in this list
		for t_seq in target_list:
			flag = True # no off-target found flag
			for i in [flank,0]:
				start_i_left = start_left+i
				w_seq_left = r["left_extend_seq"]
				if i!=0:
					w_seq_left = r["left_extend_seq"][i:-i]
				left_out = target_window_match(w_seq_left, start_i_left,t_seq,max_score=max_score)
				# out=[start_absolute,end_absolute,query_aligned,target_aligned,strand,editDistance,N_mis,N_ins,N_del,cas9_cut_position]
				start_i_right = start_right+i
				w_seq_right = r["right_extend_seq"]
				if i!=0:
					w_seq_right = r["right_extend_seq"][i:-i]
				right_out = target_window_match(w_seq_right, start_i_right,t_seq,max_score=max_score)
				# if "GACTACGTCCACCTACAGGTNGG" == t_seq:
					# print (t_seq,left_out)
					# print (t_seq,right_out)
				D_left=left_out[5]
				left_aligned=left_out[3]
				right_aligned=right_out[3]
				D_right=right_out[5]
				if D_left!=-1 and D_right!=-1:
					if D_left<=D_right:
						if left_aligned[0]=="-" or left_aligned[-1]=="-": # try longer sequence
							# print ("D_left!=-1 and D_right!=-1 1")							
							left_out2=target_window_match(r["left_extend_seq"], start_left,t_seq,max_score=max_score)
							D_left2 = left_out2[5]
							left_aligned2=left_out2[3]
							if left_aligned2[0]!="-" and left_aligned2[-1]!="-" and D_left2!=-1: #found better
								flag = False
								distance_list.append(D_left2)
								off_target_list.append([t_seq,left_out2])
								break
						flag = False
						distance_list.append(D_left)
						off_target_list.append([t_seq,left_out])
						# print (t_seq,chromosome,D_left,align1_left,align2_left)
						break
					else:
						if right_aligned[0]=="-" or right_aligned[-1]=="-": # try longer sequence 
							# print ("D_left!=-1 and D_right!=-1 2")
							right_out2=target_window_match(r["right_extend_seq"], start_right,t_seq,max_score=max_score)
							D_right2 = right_out2[5]
							right_aligned2=right_out2[3]
							if right_aligned2[0]!="-" and right_aligned2[-1]!="-" and D_right2!=-1: #found better
								flag = False
								distance_list.append(D_right2)
								off_target_list.append([t_seq,right_out2])
								break
						flag = False
						distance_list.append(D_right)
						off_target_list.append([t_seq,right_out])
						# print (t_seq,chromosome,D_right,align1_right,align2_right)
						break				
				if D_left!=-1: # not align, edit distance -1 means not align
					if left_aligned[0]=="-" or left_aligned[-1]=="-": # try longer sequence 
						# print ("D_left!=-1 and D_right!=-1 3")
						left_out2=target_window_match(r["left_extend_seq"], start_left,t_seq,max_score=max_score)
						D_left2 = left_out2[5]
						left_aligned2=left_out2[3]
						if left_aligned2[0]!="-" and left_aligned2[-1]!="-" and D_left2!=-1: #found better
							flag = False
							distance_list.append(D_left2)
							off_target_list.append([t_seq,left_out2])
							break
					flag = False
					distance_list.append(D_left)
					off_target_list.append([t_seq,left_out])
					# print (t_seq,chromosome,D_left,align1_left,align2_left)
					break
				if D_right!=-1: # not align
					if right_aligned[0]=="-" or right_aligned[-1]=="-": # try longer sequence
						# print ("D_left!=-1 and D_right!=-1 4")						
						right_out2=target_window_match(r["right_extend_seq"], start_right,t_seq,max_score=max_score)
						D_right2 = right_out2[5]
						right_aligned2=right_out2[3]
						if right_aligned2[0]!="-" and right_aligned2[-1]!="-" and D_right2!=-1: #found better
							flag = False
							distance_list.append(D_right2)
							off_target_list.append([t_seq,right_out2])
							break
					flag = False
					distance_list.append(D_right)
					off_target_list.append([t_seq,right_out])
					# print (t_seq,chromosome,D_right,align1_right,align2_right)
					break				
			if flag:
				distance_list.append(100)
		min_distance = min(distance_list)
		try:
			second_min_D = min([x for x in distance_list if x != min_distance])
		except:
			second_min_D=100
		site_list.append(min_distance)
		site_list.append(second_min_D)
		site_list.append(off_target_list)
		output_list.append(site_list)   
		del site_list
		gc.collect()
	out = pd.DataFrame(output_list,columns = df.columns.tolist()+['min_distance','second_min_D','off_target_list'])		
	del output_list
	gc.collect()
	return out


#========== MAIN ===========
file = sys.argv[1]
ncore=int(sys.argv[2])
target=sys.argv[3]
# internal paramter, same in identifyOfftargetSites.py, DO NOT CHANGE. flank is the number of bases we extend the 24bp off-target further by +/- flank. Later when we merge same off-target counts, when the +/- flank sites have the same off-target sequence, their counts will be added up.
flank = 10 
max_score = 7

#========== READ ===========
if os.path.isfile(target):
	target_list = pd.read_csv(target,header=None)[0].tolist()
else:
	target_list=[target]
target_list = [x.upper() for x in target_list]
df = pd.read_csv(file,sep="\t",header=None)
df.columns = ['#chr','start','end',"nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus","control.rev.plus","control.fwd.plus","control.fwd.minus","control.rev.minus","nuclease.mispriming.plus","nuclease.mispriming.minus","control.mispriming.plus","control.mispriming.minus",'left_extend_coord','left_extend_seq','right_extend_coord','right_extend_seq']
df.left_extend_seq = df.left_extend_seq.str.upper()
df.right_extend_seq = df.right_extend_seq.str.upper()
use_cols=['#chr','start','end',"nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus","control.rev.plus","control.fwd.plus","control.fwd.minus","control.rev.minus","nuclease.mispriming.plus","nuclease.mispriming.minus","control.mispriming.plus","control.mispriming.minus"]
#========== FILTER for faster process, ~80% rows will be filter out ===========
G_cols=["nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus"]
df2 = df[df[G_cols].sum(axis=1)<1] 
df2[use_cols].to_csv(f"{file}.noGUIDEseqReads.tsv",sep="\t",index=False)
df = df[df[G_cols].sum(axis=1)>=1]
df = df[~df.left_extend_seq.str.contains("N")]
df = df[~df.right_extend_seq.str.contains("N")]
# print (df)
def string_to_coord(df):
	df[['chr_left', 'start_left', 'end_left']] = df['left_extend_coord'].str.extract(r'(\w+):(\d+)-(\d+)')
	df[['chr_right', 'start_right', 'end_right']] = df['right_extend_coord'].str.extract(r'(\w+):(\d+)-(\d+)')
	df['start_left'] = df['start_left'].astype(int)
	df['start_right'] = df['start_right'].astype(int)
	df['end_left'] = df['end_left'].astype(int)
	df['end_right'] = df['end_right'].astype(int)
	return df
df = string_to_coord(df)
df['left_extend_seq.revcomp'] = df['left_extend_seq'].swifter.progress_bar(False).apply(reverseComplement)
df['right_extend_seq.revcomp'] = df['right_extend_seq'].swifter.progress_bar(False).apply(reverseComplement)
print ("Number of bps with >=1 GUIDE-seq-2 reads",df.shape[0])

#========== df ready for target alignment ===========
nsplit=min(200,int(df.shape[0]/2)+1)
# print (df)
dfs = np.array_split(df, nsplit)
del df
gc.collect()

df_list = Parallel(n_jobs=ncore,verbose=10)(delayed(off_target_assignment)(dfs[i],target_list,flank,max_score) for i in range(nsplit))
df = pd.concat(df_list)
# print (df.off_target_list.tolist())
del df_list
gc.collect()


# # off-target assignment
# 1. use min_distance< second_min_D
# 2. always allow multiple off-targets in one region, choose the closest cas9 as true off-target, put others into "others" column
# 2. always allow multiple off-targets in raw identified


total = df.shape[0]
df[df.min_distance>=df.second_min_D][df.columns[:19]].to_csv(f"{file}.unidentified.tsv",sep="\t",index=False)
df = df[df.min_distance<df.second_min_D]
keep = df.shape[0]

print (f"Total GUIDE-seq read start positions {total}")
print (f"Total positions with identified off-target {keep}, {keep/total}")


nsplit=min(200,int(df.shape[0]/2)+1)

dfs = np.array_split(df, nsplit)
del df
gc.collect()
def get_detailed_info(r):
	r['distance_between_cas9_cut_and_most_frequent_guideseq_read_start']=""
	r['distance_between_cas9_cut_and_most_frequent_guideseq_read_start']=abs(r['end']-r['cas9_cut_position'])
	r['total_guideseq_reads'] = r[["nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus"]].sum()
	r['total_control_reads'] = r[["control.rev.plus","control.fwd.plus","control.fwd.minus","control.rev.minus"]].sum()
	r['total_guideseq_non_dsODN_reads'] = r[["nuclease.mispriming.plus","nuclease.mispriming.minus"]].sum()
	r['total_control_non_dsODN_reads'] = r[["control.mispriming.plus","control.mispriming.minus"]].sum()
	return r


def find_off_target(df):
	# # find true off-target
	out = []
	for i,r in df.iterrows():
		off_target_list = r['off_target_list']
		current_read_start = r[df.columns[:15]].tolist()
		min_distance = r['min_distance']
		if len(off_target_list)==0:
			continue
		for t_seq,align_stat in off_target_list:
			if align_stat[5]!=min_distance:
				continue
			out.append(current_read_start+align_stat+[t_seq])
	out = pd.DataFrame(out,columns=['#chr','start','end',"nuclease.rev.plus","nuclease.fwd.plus","nuclease.fwd.minus","nuclease.rev.minus","control.rev.plus","control.fwd.plus","control.fwd.minus","control.rev.minus","nuclease.mispriming.plus","nuclease.mispriming.minus","control.mispriming.plus","control.mispriming.minus"]+["start_absolute","end_absolute","query_aligned","target_aligned","strand","editDistance","N_mis","N_ins","N_del","cas9_cut_position","on_target_sequence"])	
	# print (out)
	return out.apply(get_detailed_info,axis=1)

df_list = Parallel(n_jobs=ncore,verbose=10)(delayed(find_off_target)(dfs[i]) for i in range(nsplit))


# Note: if at A position, the read count is X, and there are multiple assigned on-targets, then the read count X is evenly assigned to each of them.
df2 = pd.concat(df_list)
del df_list
gc.collect()


df2.to_csv(f"{file}.identified.tsv",sep="\t",index=False)







