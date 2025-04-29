from __future__ import print_function
import svgwrite
import sys
import os
import logging
import argparse
import pandas as pd
import subprocess
import swifter
import numpy as np
from joblib import Parallel, delayed

logger = logging.getLogger('root')
logger.propagate = False

boxWidth = 10
box_size = 15
v_spacing = 3

# colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', 'R': '#B3B3B3', '-': '#FFFFFF'}
colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', 'R': '#B3B3B3', '-': '#B3B3B3'}
for c in ['Y','S','W','K','M','B','D','H','V','.']:
	colors[c] = "#B3B3B3"
def refseqID_to_HGNC_symbol(x,myDict):
	if "(" in x:
		ID = x.split("(")[1].split(",")[0].replace(")","").replace("(","")
		# print (ID)
		if ID in myDict:
			gene = myDict[ID]
			# print (ID,gene)
			return x.replace(ID,gene)
	return x

def reformat_homer_annotation(r):
	if r.Annotation =="Intergenic":
		return "%s (%s)"%(r.Annotation,r['Gene Name'])
	return r.Annotation
def parse_homer(identified,homer_output,genome,refseq_names=None):
	# print (refseq_names)
	out = identified.replace(".tsv",".annot.tsv")
	try:
		df = pd.read_csv(identified,sep="\t")
	except:
		return out
	select_col="Annotation"
	tmp = df[['#chr','start','end','site_name']]
	outFile = f"{identified}.homer.input"
	tmp.to_csv(outFile,sep="\t",header=False,index=False)
	command = "annotatePeaks.pl %s %s > %s"%(outFile,genome,homer_output)
	os.system(command)

	df.index = df['site_name'].to_list()
	
	df2 = pd.read_csv(homer_output,sep="\t",index_col=0)
	df2[select_col] = df2.swifter.progress_bar(False).apply(reformat_homer_annotation,axis=1)
	df['Annotation'] = df2[select_col]
	df['Annotation'] = df['Annotation'].fillna("NA")
	
	if refseq_names!=None:
		myDict = parse_HGNC(refseq_names)
		df['Annotation'] = [refseqID_to_HGNC_symbol(x,myDict) for x in df.Annotation]
	
	df.to_csv(out,sep="\t",index=False)
	return out

def get_int(x):
	try:
		x = float(x)
	except:
		return ""
	return int(x)

def parse_HGNC(f):
	refseq = "#name"
	symbol = "name2"
	df = pd.read_csv(f,sep="\t")
	# print (df.head())
	df = df[[refseq,symbol]]
	df = df.dropna()
	df.index = df[refseq].to_list()
	# print (df.head())
	return df[symbol].to_dict()
def get_alignment_info(r):
	if len(r.target_aligned)==len(r.on_target_sequence):
		r['seq'] = r.target_aligned
		r['bulged_seq'] = ""
		r['target_seq'] = r.on_target_sequence
		r['realigned_target_seq'] = ""
	elif len(r.target_aligned)>len(r.on_target_sequence): # bulge
		r['seq'] = ""
		r['bulged_seq'] = r.target_aligned
		r['target_seq'] = r.on_target_sequence
		r['realigned_target_seq'] = r.query_aligned
	else: # deletion
		r['seq'] = r.target_aligned
		r['bulged_seq'] = ""
		r['target_seq'] = r.query_aligned
		r['realigned_target_seq'] = ""
	return r
def parseSitesFile(infile):
	"""rewrite to use pandas df"""
	df = pd.read_csv(infile,sep="\t")
	df = df.sort_values("total_guideseq_reads",ascending=False)
	df['coord'] = df['#chr']+":"+df.start.astype(str)+"-"+df.end.astype(str)
	df = df.swifter.progress_bar(False).apply(get_alignment_info,axis=1)
	# df['seq'] = df['target_aligned']
	# df['bulged_seq'] = df['Site_GapsAllowed.site_sequence'].fillna("")
	# df['target_seq'] = df['on_target_sequence']
	# df['realigned_target_seq'] = df['Site_GapsAllowed.on_target_sequence'].replace("none","")
	df['reads'] = df['total_guideseq_reads'].astype(str)
	df['num_mismatch'] = df['editDistance'].astype(int).astype(str)
	try:
		df['annot'] = df['Annotation']
	except:
		df['annot'] = ""
	try:
		df['Title'] = df['Title']
	except:
		df['Title'] = ""
	try:
		df['outfile'] = df['outfile']
	except:
		df['outfile'] = ""

	outDict=df.groupby('on_target_sequence')[['seq', 'bulged_seq', 'reads', 'coord', 'annot', 'num_mismatch', 'target_seq', 'realigned_target_seq','Title','outfile']].apply(lambda x: x.to_dict(orient='records')).to_dict()

	
	return outDict
def parseSitesFile_bk(infile):
	"""rewrite to use pandas df"""
	df = pd.read_csv(infile,sep="\t")
	df = df.sort_values("total_guideseq_reads",ascending=False)
	df['coord'] = df['#site_chr']+":"+df.site_start.astype(str)+"-"+df.site_end.astype(str)
	df['seq'] = df['Site_SubstitutionsOnly.site_sequence'].fillna("")
	df['bulged_seq'] = df['Site_GapsAllowed.site_sequence'].fillna("")
	df['reads'] = df['total_guideseq_reads'].astype(str)
	df['num_mismatch'] = np.where(df['Site_SubstitutionsOnly.distance'] != "", 
			df['Site_SubstitutionsOnly.distance'], 
			df['Site_GapsAllowed.distance'])
	df['num_mismatch'] = df['num_mismatch'].astype(int).astype(str)
	try:
		df['annot'] = df['Annotation']
	except:
		df['annot'] = ""
	try:
		df['Title'] = df['Title']
	except:
		df['Title'] = ""
	try:
		df['outfile'] = df['outfile']
	except:
		df['outfile'] = ""
	df['target_seq'] = df['on_target_sequence']
	df['realigned_target_seq'] = df['Site_GapsAllowed.on_target_sequence'].replace("none","")
	outDict=df.groupby('on_target_sequence')[['seq', 'bulged_seq', 'reads', 'coord', 'annot', 'num_mismatch', 'target_seq', 'realigned_target_seq','Title','outfile']].apply(lambda x: x.to_dict(orient='records')).to_dict()

	
	return outDict

# 3/6/2020

from Bio import SeqUtils
def find_PAM(seq,PAM):
	try:
		PAM_index = seq.index(PAM)
							 
							 
															
	except:
		# PAM on the left
		left_search = SeqUtils.nt_search(seq[:len(PAM)], PAM)
		if len(left_search)>1:
			PAM_index = left_search[1]
		else:
			right_search = SeqUtils.nt_search(seq[-len(PAM):], PAM)
			if len(right_search)>1:
				PAM_index = len(seq)-len(PAM)
			else:
				print ("PAM: %s not found in %s. Set PAM index to 20"%(PAM,seq))
				PAM_index=20
	return PAM_index

def draw_svg(target_seq,PAM,offtargets,title,genome,outfile):
	# Initiate canvas
	outdir = os.path.dirname(outfile)
	# print (outdir)
	if offtargets[0]['outfile']!="":
		dwg = svgwrite.Drawing(f"{outdir}/{offtargets[0]['outfile']}.svg", profile='full', size=(u'100%', 100 + (len(offtargets)+5)*(box_size + 1)))
	else:
		dwg = svgwrite.Drawing(f"{outfile}.svg", profile='full', size=(u'100%', 100 + (len(offtargets)+5)*(box_size + 1)))
	if offtargets[0]['Title']!="":
		# Define top and left margins
		x_offset = 20
		y_offset = 50
		dwg.add(dwg.text(offtargets[0]['Title'], insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
	elif title is not None:
		# Define top and left margins
		x_offset = 20
		y_offset = 50
		dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
	else:
		# Define top and left margins
		x_offset = 20
		y_offset = 20

	## Assume PAM is on the right end Yichao rewrite visualization code, generic PAM
	## PAM can be on the left or right, Yichao 0713
	tick_locations = []
	tick_legend = []
	# PAM_index = target_seq.index(PAM)
	PAM_index = find_PAM(target_seq,PAM)
	count = 0
	for i in range(PAM_index,0,-1):
		count = count+1
		if count % 10 == 0:
			tick_legend.append(count)
			# print (count,i)
			tick_locations.append(i)
	if len(PAM)>=3:
		tick_legend+=['P', 'A', 'M']+['-']*(len(PAM)-3)
	else:
		tick_legend+=["PAM"]+['-']*(len(PAM)-3)
	tick_locations+=range(PAM_index+1,len(target_seq)+1)
	if PAM_index == 0:
		tick_legend = []
		tick_locations = []
		tick_legend+=['P', 'A', 'M']+['-']*(len(PAM)-3)
		tick_locations+=range(1,len(PAM)+1)
		count = 0
		for i in range(len(PAM)+1,len(target_seq)+1):
			count = count+1
			if count % 10 == 0 or count == 1:
				tick_legend.append(count)
				# print (count,i)
				tick_locations.append(i)
	# print (zip(tick_locations, tick_legend))
	for x,y in zip(tick_locations, tick_legend):
		dwg.add(dwg.text(y, insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

	# Draw reference sequence row
	for i, c in enumerate(target_seq):
		y = y_offset
		x = x_offset + i * box_size
		dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
		dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))
	# dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(target_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
	# dwg.add(dwg.text('Mismatches', insert=(box_size * (len(target_seq) + 1) + 90, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
	# dwg.add(dwg.text('Coordinates', insert=(box_size * (len(target_seq) + 1) + 200, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
	# if genome!=None:
		# dwg.add(dwg.text('Annotation', insert=(box_size * (len(target_seq) + 1) + 450, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
	dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(target_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
	dwg.add(dwg.text('Mismatches', insert=(box_size * (len(target_seq) + 1) + 150, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
	dwg.add(dwg.text('Coordinates', insert=(box_size * (len(target_seq) + 1) + 250, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
	if genome!=None:
		dwg.add(dwg.text('Annotation', insert=(box_size * (len(target_seq) + 1) + 500, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))

	# Draw aligned sequence rows
	y_offset += 1  # leave some extra space after the reference row
	line_number = 0  # keep track of plotted sequences

										
	for j, seq in enumerate(offtargets):
		realigned_target_seq = offtargets[j]['realigned_target_seq']
		no_bulge_offtarget_sequence = offtargets[j]['seq']
		bulge_offtarget_sequence = offtargets[j]['bulged_seq']

		if no_bulge_offtarget_sequence != '':
			k = 0
			line_number += 1
			y = y_offset + line_number * box_size
			for i, (c, r) in enumerate(zip(no_bulge_offtarget_sequence, target_seq)):
				x = x_offset + k * box_size
				if r == '-':
					if 0 < k < len(target_seq):
						x = x_offset + (k - 0.25) * box_size
						dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
						dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
				elif c == r:
					dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
					k += 1
				elif r == 'N':
					dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
					k += 1
				else:
					dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
					dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
					k += 1
		if bulge_offtarget_sequence != '':
			k = 0
			line_number += 1
			y = y_offset + line_number * box_size
			for i, (c, r) in enumerate(zip(bulge_offtarget_sequence, realigned_target_seq)):
				x = x_offset + k * box_size
				if r == '-':
					if 0 < k < len(realigned_target_seq):
						x = x_offset + (k - 0.25) * box_size
						dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
						dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
				elif c == r:
					dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
					k += 1
				elif r == 'N':
					dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
					k += 1
				else:
					dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
					dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
					k += 1

		if no_bulge_offtarget_sequence == '' or bulge_offtarget_sequence == '':
			reads_text = dwg.text(str(seq['reads']), insert=(box_size * (len(target_seq) + 1) + 20, y_offset + box_size * (line_number + 2) - 2),
								  fill='black', style="font-size:15px; font-family:Courier")
			dwg.add(reads_text)
			mismatch_text = dwg.text(seq['num_mismatch'], insert=(box_size * (len(target_seq) + 1) + 180, y_offset + box_size * (line_number + 2) - 2),
								  fill='black', style="font-size:15px; font-family:Courier")
			dwg.add(mismatch_text)
			mismatch_text = dwg.text(seq['coord'], insert=(box_size * (len(target_seq) + 1) + 250, y_offset + box_size * (line_number + 2) - 2),
								  fill='black', style="font-size:15px; font-family:Courier")
			dwg.add(mismatch_text)
			if genome!= None:
				annot_text = dwg.text(seq['annot'], insert=(box_size * (len(target_seq) + 1) + 500, y_offset + box_size * (line_number + 2) - 2),
									  fill='black', style="font-size:15px; font-family:Courier")
				dwg.add(annot_text)
		else:
			reads_text = dwg.text(str(seq['reads']), insert=(box_size * (len(target_seq) + 1) + 20, y_offset + box_size * (line_number + 1) + 5),
								  fill='black', style="font-size:15px; font-family:Courier")
			dwg.add(reads_text)
			mismatch_text = dwg.text(seq['num_mismatch'], insert=(box_size * (len(target_seq) + 1) + 180, y_offset + box_size * (line_number + 1) + 5),
								  fill='black', style="font-size:15px; font-family:Courier")
			dwg.add(mismatch_text)
			mismatch_text = dwg.text(seq['coord'], insert=(box_size * (len(target_seq) + 1) + 250, y_offset + box_size * (line_number + 1) + 5),
								  fill='black', style="font-size:15px; font-family:Courier")
			dwg.add(mismatch_text)
			if genome!= None:
				annot_text = dwg.text(seq['annot'], insert=(box_size * (len(target_seq) + 1) + 500, y_offset + box_size * (line_number + 1) + 5),
									  fill='black', style="font-size:15px; font-family:Courier")
				dwg.add(annot_text)
			reads_text02 = dwg.text(u"\u007D", insert=(box_size * (len(target_seq) + 1) + 7, y_offset + box_size * (line_number + 1) + 5),
								  fill='black', style="font-size:23px; font-family:Courier")
			dwg.add(reads_text02)
	dwg.save()



def visualizeOfftargets(infile, outfile, title, PAM, genome=None,refseq_names=None,njobs=1):

	output_folder = os.path.dirname(outfile)
	if not os.path.exists(output_folder):
		try:
			os.makedirs(output_folder)
		except:
			print ("Folder exist")

	if genome!=None:
		infile = parse_homer(infile,outfile+".raw.homer.tsv",genome,refseq_names=refseq_names)

	myDict = parseSitesFile(infile)
	Parallel(n_jobs=njobs,verbose=10)(delayed(draw_svg)(target_seq,PAM,myDict[target_seq],title,genome,f"{outfile}_{target_seq}") for target_seq in myDict)

	# for target_seq in myDict:
		# print (target_seq,myDict[target_seq])
		# draw_svg(target_seq,PAM,myDict[target_seq],title,genome,f"{outfile}_{target_seq}")
	return myDict
	
def visualizeOfftargets_skip_homer(infile, outfile, title, PAM, genome=None,refseq_names=None,njobs=1):

	output_folder = outfile
	if not os.path.exists(output_folder):
		try:
			os.makedirs(output_folder)
		except:
			print ("Folder exist")

	myDict = parseSitesFile(infile)
	Parallel(n_jobs=njobs,verbose=10)(delayed(draw_svg)(target_seq,PAM,myDict[target_seq],title,genome,f"{outfile}/{outfile}_{target_seq}") for target_seq in myDict)

	# for target_seq in myDict:
		# print (target_seq,myDict[target_seq])
		# draw_svg(target_seq,PAM,myDict[target_seq],title,genome,f"{outfile}_{target_seq}")
	return myDict

def main():
	parser = argparse.ArgumentParser(description='Plot visualization plots for re-aligned reads.')
	parser.add_argument("-f","--identified_file", help="", required=True)
	parser.add_argument("-o","--outfile", help="FullPath/VIZ", required=True)
	parser.add_argument("-t","--title", help="Plot title", required=True)
	parser.add_argument("--skip_homer",  help="skip homer",action='store_true')		
	parser.add_argument("-g","--genome", help="if specified, homer annotation will be performed", default=None)
	parser.add_argument("-n","--njobs", help="number of cores", default=1,type=int)
	parser.add_argument("-a","--annotation", help="refseqID to gene name mapping", default=None)
	parser.add_argument("--PAM", help="PAM sequence", default="NGG")	
	args = parser.parse_args()

	print(args)
	if args.skip_homer:
		visualizeOfftargets_skip_homer(args.identified_file, args.outfile, args.title, args.PAM,args.genome,args.annotation,args.njobs)
	else:
		visualizeOfftargets(args.identified_file, args.outfile, args.title, args.PAM,args.genome,args.annotation,args.njobs)

if __name__ == "__main__":

	main()