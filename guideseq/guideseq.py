#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

guideseq.py
===========

guideseq pool

serves as the wrapper for all guideseq pipeline

"""
import glob
import pandas as pd
import os
import sys
import yaml
import argparse
import traceback
import subprocess
# Set up logger
import log
logger = log.createCustomLogger('root',"guideseq_pool")
import datetime
from alignReads import alignReads
from filterBackgroundSites import filterBackgroundSites,filterBlackList,filterControl
from visualization import visualizeOfftargets
import identifyOfftargetSites
import validation
from tabulate import tabulate
from demultiplex import demultiplex,demultiplex_parallel
from joblib import Parallel, delayed

DEFAULT_YAML = os.path.dirname(os.path.realpath(__file__)) + "/default.yaml"
myDate=str(datetime.date.today())
def get_parameters(manifest_data):

	# init
	default_refseqName = os.path.dirname(os.path.realpath(__file__)) + "/refseq_gene_name.py"
	with open(DEFAULT_YAML, 'r') as f:
		default = yaml.load(f)
	with open(manifest_data, 'r') as f:
		return_dict = yaml.load(f) # this is user input YAML
	default['analysis_folder'] = os.getcwd()
	default['refseq_names'] = default_refseqName
	validation.validateManifest(return_dict)
	
	# assign default if user YAML missing some parameters
	return_dict = assign_default(return_dict,default)

	return return_dict

def assign_default(return_dict,default):
	
	section_list = ['demultiplex']
	for k in section_list:
		if k in return_dict:
			for p in default[k]:
				if not p in return_dict[k]:
					return_dict[k][p] = default[k][p]
		else:
			return_dict[k] = default[k]

	# default global variables
	for p in default:
		if not p in section_list+['samples']:
			if not p in return_dict:
				return_dict[p] = default[p]
	
	# set abs path for input and output
	for k in section_list:
		return_dict[k]['out_dir'] = os.path.join(return_dict['analysis_folder'], return_dict[k]['out_dir'])
		return_dict[k]['input_dir'] = os.path.join(return_dict['analysis_folder'], return_dict[k]['input_dir'])

	return return_dict

class GuideSeq:

	def __init__(self):
		pass

	def parseManifest(self, manifest_path,sample='all',user_dict=None):
		logger.info('Loading manifest...')

		try:
			manifest_data = get_parameters(manifest_path)
			self.samples = {}
			self.parameters = manifest_data
			# print (self.parameters['guideseq_pool']==False)
			# exit()
			if user_dict == "None":
				user_dict = None
			if user_dict:
				for p in user_dict:
					self.parameters[p] = user_dict[p]
			if sample != "all":
				self.samples[sample] = manifest_data['samples'][sample]
			else:
				self.samples = manifest_data['samples']
			del self.parameters['samples']
			logger.info("\n"+tabulate([(k,v) for k,v in self.parameters.items()])) 
			logger.info("\n"+tabulate([(k,v) for k,v in self.samples[list(self.samples.keys())[0]].items()])) 

		except Exception as e:
			logger.error(
				'Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
			logger.error(traceback.format_exc())
			sys.exit()

	def demultiplex(self):

		barcode_dict = {}
		for sample in self.samples:
			barcode1 = self.samples[sample]['barcode1']
			barcode2 = self.samples[sample]['barcode2']
			barcode_dict[sample] = [barcode1,barcode2]
			
			barcode1 = self.samples[sample]['controlbarcode1']
			barcode2 = self.samples[sample]['controlbarcode2']
			barcode_dict["Control_"+sample] = [barcode1,barcode2]
		# print (barcode_dict)
		try:
			logger.info('Demultiplexing Undetermined files...')
			if self.parameters['njobs'] == 1 or self.parameters['guideseq_pool']==False:
				logger.info('guide-seq demultiplex mode, 1-core')
				demultiplex(os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['forward']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['reverse']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['index1']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['index2']),
						barcode_dict,
						out_dir=self.parameters['demultiplex']['out_dir'],
						mismatch=self.parameters['demultiplex']['mismatch'],
						min_reads=self.parameters['demultiplex']['min_reads'],subsample_reads=self.parameters['demultiplex']['subsample_reads'],revcomp_I2=self.parameters['demultiplex']['revcomp_I2'])
			else:
				demultiplex_parallel(os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['forward']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['reverse']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['index1']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['index2']),
						barcode_dict,
						out_dir=self.parameters['demultiplex']['out_dir'],
						mismatch=self.parameters['demultiplex']['mismatch'],
						min_reads=self.parameters['demultiplex']['min_reads'],subsample_reads=self.parameters['demultiplex']['subsample_reads'],splitFastq_path=self.parameters['splitFastq'],ncore=self.parameters['njobs'],revcomp_I2=self.parameters['demultiplex']['revcomp_I2'])

			logger.info('Successfully demultiplexed reads.')

		except Exception as e:
			logger.error('Error demultiplexing reads.')
			logger.error(traceback.format_exc())
			quit()

	def alignReads(self):
		logger.info('Aligning reads and collapsing reads based on UMI...')

		try:
			for sample in self.samples:
				try:
					# default input
					if not "read1" in self.samples[sample]:
						if self.parameters['demultiplex']['subsample_reads'] > 0:
							self.samples[sample]['read1'] = os.path.join(self.parameters['demultiplex']['out_dir'], f"{sample}.sampling.r1.fastq")
							self.samples[sample]['read2'] = os.path.join(self.parameters['demultiplex']['out_dir'], f"{sample}.sampling.r2.fastq")
							self.samples[sample]['controlread1'] = os.path.join(self.parameters['demultiplex']['out_dir'], f"Control_{sample}.sampling.r1.fastq")
							self.samples[sample]['controlread2'] = os.path.join(self.parameters['demultiplex']['out_dir'], f"Control_{sample}.sampling.r2.fastq")
						else:
							self.samples[sample]['read1'] = os.path.join(self.parameters['demultiplex']['out_dir'], f"{sample}.r1.fastq")
							self.samples[sample]['read2'] = os.path.join(self.parameters['demultiplex']['out_dir'], f"{sample}.r2.fastq")
							self.samples[sample]['controlread1'] = os.path.join(self.parameters['demultiplex']['out_dir'], f"Control_{sample}.r1.fastq")
							self.samples[sample]['controlread2'] = os.path.join(self.parameters['demultiplex']['out_dir'], f"Control_{sample}.r2.fastq")
					
					sample_alignment_path = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.dedup.bam')
					# def alignReads(bwa=None, samtools=None, # Tools
					# reference_genome=None, R1=None, R2=None, label=None, output_dir=None, # Inputs
					# ,umi_tools=None, njobs=None, **kwargs)
					alignReads(
							   R1=self.samples[sample]['read1'],
							   R2=self.samples[sample]['read2'],
							   label=sample,
							   output_dir=f"{self.parameters['analysis_folder']}/aligned",
							   **self.parameters
							   )
					# exit()
					self.samples[sample]['aligned'] = sample_alignment_path
					sample_alignment_path = os.path.join(self.parameters['analysis_folder'], 'aligned', "Control_"+sample + '.dedup.bam')
					alignReads(
							   R1=self.samples[sample]['controlread1'],
							   R2=self.samples[sample]['controlread2'],
							   label="Control_"+sample,
							   output_dir=f"{self.parameters['analysis_folder']}/aligned",
							   **self.parameters
							   )
					self.samples[sample]['controlaligned'] = sample_alignment_path
				except:
					logger.error(f'Failed for sample {sample}')	
					logger.error(traceback.format_exc())
				logger.info('Finished aligning reads to genome.')

		except Exception as e:
			logger.error('Error aligning')
			logger.error(traceback.format_exc())
			quit()

	def identifyOfftargetSites(self):
		logger.info('Identifying offtarget sites...')

		try:
			self.identified = {}

			# Identify offtarget sites for each sample
			for sample in self.samples:
				try:
					# default input
					if not "aligned" in self.samples[sample]:
						try:
							self.samples[sample]['aligned'] = self.parameters['aligned']
							self.samples[sample]['controlaligned'] = self.parameters['controlaligned']
							logger.info('Using user provided dedup BAM path...')
						except:
							logger.info('Using system default dedup BAM path...')
							self.samples[sample]['aligned'] = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.dedup.bam')
							self.samples[sample]['controlaligned'] = os.path.join(self.parameters['analysis_folder'], 'aligned', "Control_"+sample + '.dedup.bam')
					# Prepare sample annotations
					sample_data = self.samples[sample]

					self.identified[sample] = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '.matched.final.high_confidence.tsv')
					# new function
					identifyOfftargetSites.identify_sites(self.parameters['reference_genome'],
					   self.samples[sample]['aligned'],
					   self.samples[sample]['controlaligned'],
					   self.samples[sample]['aligned'].replace(f"aligned/{sample}",f"aligned/{sample}.mispriming"),
					   self.samples[sample]['controlaligned'].replace(f"aligned/Control_{sample}",f"aligned/Control_{sample}.mispriming"),
					   sample,
					   f"{self.parameters['analysis_folder']}/identified",
					   self.samples[sample]['target'],
					   njobs=self.parameters['njobs'],
					   extend_spacer=self.parameters['extend_spacer'],
					   extend_pam=self.parameters['extend_pam'],
					   flank=self.parameters['extend_flank'],
					   bedtools=self.parameters['bedtools'],
					   chrom_size=self.parameters['chrom_size'],
					   min_read_count=self.parameters['min_guideseq_read_count'], # not used
					   max_control_reads=self.parameters['max_control_read_count'] # not used
					   )

				except:
					logger.error(f'Failed for sample {sample}')   
					logger.error(traceback.format_exc())

			logger.info('Finished identifying offtarget sites.')

		except Exception as e:
			logger.error('Error identifying offtarget sites.')
			logger.error(traceback.format_exc())
			quit()

	def visualize(self):
		logger.info('Visualizing off-target sites')


		for sample in self.samples: ## 3/6/2020 Yichao solved: visualization stopped when one sample failed
			try:
				infile = self.identified[sample]
				outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '_offtargets')
				logger.info("visualizeOfftargets")
				visualizeOfftargets(infile, outfile, title=sample,PAM=self.parameters['PAM'],genome=self.parameters['genome'],refseq_names=self.parameters['refseq_names'],njobs=self.parameters['njobs'])

			except Exception as e:
				logger.error('Error visualizing off-target sites: %s'%(sample))
				logger.error(traceback.format_exc())

		logger.info('Finished visualizing off-target sites')
	def parallel(self, manifest_path, lsf,step,overwrite,queue):
		logger.info('Submitting parallel jobs for GUIDEseq_pool')
		if not os.path.exists("HPC_parallel_log"):
			os.makedirs("HPC_parallel_log")
		current_script = __file__
		count = 1
		if "gpu" in queue:
			if overwrite == None:
				overwrite={}
			overwrite['bwa']='fq2bam'
			lsf+= ' -gpu "num=1/host:mode=exclusive_process"'
		try:
			for sample in self.samples:
				cmd = f'/research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/anaconda3/envs/changeseq_py3/bin/python {current_script} main --manifest {manifest_path} --sample {sample} --step {step} --overwrite "{overwrite}"'
				cmd=lsf.replace("{Sample_Name}",sample).replace("{queue}",queue).replace("{myDate}",myDate).replace("{memory}",str(self.parameters['memory'])).replace("{njobs}",str(self.parameters['njobs'])) + f" -J {sample[:10]} " + cmd
				logger.info(cmd)
				subprocess.call(cmd,shell=True)
				count += 1
			logger.info('Finished job submission')

		except Exception as e:
			logger.error('Error submitting jobs.')
			logger.error(traceback.format_exc())


def parse_args():
	parser = argparse.ArgumentParser()

	subparsers = parser.add_subparsers(description='Two run mode: main or parallel',
									   help='Use this to run individual steps of the pipeline',
									   dest='command')

	parallel_parser = subparsers.add_parser('parallel', help='submit a job for each sample specified in the YAML file. Tested in the LSF system. May not work in other job management systems.')
	parallel_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	parallel_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')
	parallel_parser.add_argument('--queue', '-q', help='LSF queue', default='standard')
	parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -R "rusage[mem={memory}] span[hosts=1]" -n {njobs} -P GUIDEseq_pool -q {queue} -o HPC_parallel_log/GUIDEseq_pool_{myDate}_{Sample_Name}_%J.log')
	parallel_parser.add_argument('--step', help='Specify which steps of pipepline to run (demultiplex, align, identify,visualize)', default='demultiplex+align+identify+visualize')
	# parallel_parser.add_argument('--step', help='Specify which steps of pipepline to run (demultiplex, align, identify,visualize)', default='demultiplex+align')
	parallel_parser.add_argument('--overwrite', help='overwrite specifications in the yaml file', default=None,type=yaml.load)

	main_parser = subparsers.add_parser('main', help='Run a single step or a series of steps for one or all samples')
	main_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	main_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')
	main_parser.add_argument('--step', help='Specify steps, demultiplex, align, identify,visualize, order does not matter', default='demultiplex+align+identify+visualize')
	main_parser.add_argument('--overwrite', help='overwrite specifications in the yaml file', default=None,type=yaml.load)
	
	QC_parser = subparsers.add_parser('qc', help='summerize QC')
	QC_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)

	
	init_parser = subparsers.add_parser('init', help='initialize a default yaml in the current folder')

	return parser.parse_args()


def main():
	args = parse_args()
	

	if args.command == "init":
		subprocess.call(f"cp {DEFAULT_YAML} Guideseq.yaml",shell=True)
		print("Guideseq.yaml has been initialize. \nPlease add sample info and then run GUIDE-seq analysis using the 'main' or the 'parallel' subcommand.")
	elif str(args.command).lower() == "qc":
		logger.info("Running QC stats")
		g = GuideSeq()
		g.parseManifest(args.manifest)
		output_dir=g.parameters['analysis_folder']+"/QC"
		if not os.path.exists(output_dir):
			try:
				os.makedirs(output_dir)
				print ("creating folder",output_dir)
			except:
				print ("Folder exist")
		
		command_list = []
		# fastqc
		command_list.append(f"fastqc {os.path.join(g.parameters['demultiplex']['input_dir'],g.parameters['demultiplex']['forward'])} -o {output_dir}")
		command_list.append(f"fastqc {os.path.join(g.parameters['demultiplex']['input_dir'],g.parameters['demultiplex']['reverse'])} -o {output_dir}")
		command_list.append(f"fastqc {os.path.join(g.parameters['demultiplex']['input_dir'],g.parameters['demultiplex']['index1'])} -o {output_dir}")
		command_list.append(f"fastqc {os.path.join(g.parameters['demultiplex']['input_dir'],g.parameters['demultiplex']['index2'])} -o {output_dir}")
		files = glob.glob(f"{g.parameters['analysis_folder']}/*/*fastq")
		for f in files:
			command_list.append(f'fastqc {f} -o {output_dir}')
		
		# flagstat
		files = glob.glob(f"{g.parameters['analysis_folder']}/aligned/*st.bam")+glob.glob(f"{g.parameters['analysis_folder']}/aligned/*dedup.bam")
		for f in files:
			command_list.append(f'samtools flagstat {f} > {f}.flagstat')
		
		
		Parallel(n_jobs=g.parameters['njobs'],verbose=10)(delayed(os.system)(command) for command in command_list)
		#
		out = []
		for sample in g.samples:
			f = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '.matched.final.high_confidence.tsv')
			if os.path.isfile(f):
				df = pd.read_csv(f,sep="\t")
				for s, d in df.groupby("on_target_sequence"):
					on_target = d[d.editDistance==0]
					total_reads = d.total_guideseq_reads.sum()
					on_target_reads = on_target.total_guideseq_reads.sum()
					N_on_target = on_target.shape[0]
					N_off_target = d.shape[0]-N_on_target
					out.append([sample,s,on_target_reads,N_on_target,N_off_target])
		out = pd.DataFrame(out,columns=['Sample','Target','#Reads_on_target','#N_on_target','#N_off_target'])
		out.Sample = out.Sample+"_"+out.Target
		out = out.drop(['Target'],axis=1)
		outfile=f"{os.path.join(g.parameters['analysis_folder'], 'identified', 'identified.stats_mqc.tsv')}"
		with open(outfile, 'w') as f:
			f.write("# plot_type: 'table'\n")
			f.write("# section_name: 'GUIDE-seq on-target STATS'\n")
		out.to_csv(outfile,
				   sep="\t",
				   index=False,
				   mode='a',
				   header=True)  # or header=True if you want column names after the comments

		
	elif args.command == "main":
		# run a single sample
		logger.info("User command: %s"%(" ".join(sys.argv)))
		g = GuideSeq()
		g.parseManifest(args.manifest,args.sample,args.overwrite)
		steps  = args.step.split("+")
		if "demultiplex" in steps:
			g.demultiplex()
		if "align" in steps:
			g.alignReads()
		if "identify" in steps:
			g.identifyOfftargetSites()
			# exit()
			# g.filterBackgroundSites()
		if "visualize" in steps:
			try:
				g.identified
			except:
				g.identified = {}
				for sample in g.samples:
					if 'identified' in g.samples[sample]:
						g.identified[sample] = g.samples[sample]['identified']
					else:
						g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '.matched.final.high_confidence.tsv')
			g.visualize()

	elif args.command == 'parallel':
		c = GuideSeq()
		c.parseManifest(args.manifest,args.sample,args.overwrite)
		# exit()
		c.parallel(args.manifest, args.lsf, args.step,args.overwrite,args.queue)

 
if __name__ == '__main__':
	main()
