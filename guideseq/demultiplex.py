from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import time
import logging
from Levenshtein import distance,hamming
import subprocess
from joblib import Parallel, delayed
import glob
logger = logging.getLogger('root')
import os
import sys
def fq(file):
	if re.search('.gz$', file):
		fastq = gzip.open(file, 'rt')
	else:
		fastq = open(file, 'r')
	with fastq as f:
		while True:
			l1 = f.readline()
			if not l1:
				break
			l2 = f.readline()
			l3 = f.readline()
			l4 = f.readline()
			yield [l1, l2, l3, l4]


def get_sample_id(i1, i2, sample_names,mismatch=1):
	seq1 = i1[1]
	seq2 = i2[1]
	# print (seq1[:8],seq2[:8])
	for k in sample_names:
		b1,b2 = sample_names[k]
		if distance(seq1[:8],b1) <= mismatch:
			if distance(seq2[:8],b2) <= mismatch:
				return k
	return "skip"

def demultiplex_parallel(read1, read2, index1, index2, sample_names={}, out_dir="out_dir",mismatch=1,min_reads=0,subsample_reads=-1,splitFastq_path=None,ncore=None):
	if os.system(splitFastq_path) == 0:
		print ("splitFastq needs to be installed to run demultiplex in parallel, running in single core mode")
		demultiplex(read1, read2, index1, index2, sample_names, out_dir,mismatch,min_reads,subsample_reads)
		return 1
	if not os.path.exists(out_dir):
		try:
			os.makedirs(out_dir)
			print ("creating folder",out_dir)
		except:
			print ("Folder exist")
	number_reads=1000000
	command = f"{splitFastq_path} -n {number_reads} -i {read1} -o {out_dir}/Undetermined_split.R1"
	subprocess.call(command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
	command = f"{splitFastq_path} -n {number_reads} -i {read2} -o {out_dir}/Undetermined_split.R2"
	subprocess.call(command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
	command = f"{splitFastq_path} -n {number_reads} -i {index1} -o {out_dir}/Undetermined_split.I1"
	subprocess.call(command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
	command = f"{splitFastq_path} -n {number_reads} -i {index2} -o {out_dir}/Undetermined_split.I2"
	subprocess.call(command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
	R1files = glob.glob(f"{out_dir}/*Undetermined_split.R1.fastq")
	# print (R1files)
	Parallel(n_jobs=ncore,verbose=10)(delayed(demultiplex)(R1, R1.replace("R1","R2"), R1.replace("R1","I1"), R1.replace("R1","I2"), sample_names, R1.split(".R1.")[0]+".split_demultiplex_out",mismatch,0,-1) for R1 in R1files)
	# merge
	for sample_id in sample_names:
		command = f"cat {out_dir}/*split_demultiplex_out/{sample_id}*.r1.fastq > {out_dir}/{sample_id}.r1.fastq"
		subprocess.call(command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
		command = f"cat {out_dir}/*split_demultiplex_out/{sample_id}*.r2.fastq > {out_dir}/{sample_id}.r2.fastq"
		subprocess.call(command,shell=True,stdout=sys.stdout,stderr=sys.stderr)


def demultiplex(read1, read2, index1, index2, sample_names={}, out_dir="out_dir",mismatch=1,min_reads=10000,subsample_reads=-1):
	"""demultiplexing guideseq Undetermined fastq
	
	Input
	-----
	
	sample_names['sample_id'] = [barcode1, barcode2]
	
	"""
	if not os.path.exists(out_dir):
		try:
			os.makedirs(out_dir)
			print ("creating folder",out_dir)
		except:
			print ("Folder exist")
	
	outfiles_r1 = {}
	outfiles_r2 = {}
	# outfiles_i1 = {}
	# outfiles_i2 = {}

	total_count = 0
	count = {}
	# print (read1,read2)

	start = time.time()
	for r1,r2,i1,i2 in zip(fq(read1), fq(read2), fq(index1), fq(index2)):
		total_count += 1
		if total_count % 1000000 == 0:
			logger.info("Processed %d reads in %.1f minutes.", total_count, (time.time()-start)/60)
		sample_id = get_sample_id(i1, i2, sample_names,mismatch=mismatch)
		if sample_id == "skip":
			continue
		# Increment read count and create output buffers if this is a new sample barcode
		if not sample_id in count:
			count[sample_id] = 0
		count[sample_id] += 1

		if not sample_id in outfiles_r1:
			outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % sample_id), 'w')
			outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % sample_id), 'w')
		name,barcode = r1[0].strip().split()
		umi = i2[1].strip()[-10:]+r1[1][:5]+r2[1][:5]
		# umi = barcode.split("+")[-1][-10:]
		r1[0] = name+"_"+umi+"\n"
		name,barcode = r2[0].strip().split()
		# umi = barcode.split("+")[-1][-10:]
		r2[0] = name+"_"+umi+"\n"
		for line in r1:
			print (line, file=outfiles_r1[sample_id], end="")
		for line in r2:
			print (line, file=outfiles_r2[sample_id], end="")
			# for line in i1:
				# print (line, file=outfiles_i1[sample_id], end="")
			# for line in i2:
				# print (line, file=outfiles_i2[sample_id], end="")

	# Write remaining buffered reads to a single fastq.
	# (These reads correspond to barcodes that were seen less than min_reads times)
	# undetermined_r1 = open(os.path.join(out_dir, 'undetermined.r1.fastq'), 'w')
	# undetermined_r2 = open(os.path.join(out_dir, 'undetermined.r2.fastq'), 'w')
	# undetermined_i1 = open(os.path.join(out_dir, 'undetermined.i1.fastq'), 'w')
	# undetermined_i2 = open(os.path.join(out_dir, 'undetermined.i2.fastq'), 'w')
	# for sample_id in buffer_r1.keys():
		# for record in buffer_r1[sample_id]:
			# undetermined_r1.write(''.join(record))
		# for record in buffer_r2[sample_id]:
			# undetermined_r2.write(''.join(record))
		# for record in buffer_i1[sample_id]:
			# undetermined_i1.write(''.join(record))
		# for record in buffer_i2[sample_id]:
			# undetermined_i2.write(''.join(record))

	# Close files
	for sample_id in outfiles_r1:
		outfiles_r1[sample_id].close()
		outfiles_r2[sample_id].close()
		# outfiles_i1[sample_id].close()
		# outfiles_i2[sample_id].close()
	# undetermined_r1.close()
	# undetermined_r2.close()
	# undetermined_i1.close()
	# undetermined_i2.close()
	try:
		num_fastqs = len([v for k,v in count.items() if v>=min_reads])
		logger.info('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads.', num_fastqs, len(count), min_reads)
		min_reads = min(count.values())
		# logger.info(min_reads)
		if subsample_reads > 0:
			logger.info("subsampling reads...")
			if subsample_reads > min_reads:
				logger.info('User input subsample_reads is higher than min reads')
				logger.info('setting subsample_reads to: %s'%(min_reads))
				subsample_reads = min_reads-1
			for sample_id in count:
				input_R1 = os.path.join(out_dir, '%s.r1.fastq' % sample_id)
				input_R2 = os.path.join(out_dir, '%s.r2.fastq' % sample_id)
				out_R1 = os.path.join(out_dir, '%s.sampling.r1.fastq' % sample_id)
				out_R2 = os.path.join(out_dir, '%s.sampling.r2.fastq' % sample_id)
				command_R1 = "seqtk sample -s100 %s %s > %s"%(input_R1,subsample_reads,out_R1)
				command_R2 = "seqtk sample -s100 %s %s > %s"%(input_R2,subsample_reads,out_R2)
				os.system(command_R1)
				os.system(command_R2)
	except Exception as e:
		print (e)
		print (count)
