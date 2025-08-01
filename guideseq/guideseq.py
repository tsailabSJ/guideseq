#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

guideseq.py
===========

V2

serves as the wrapper for all guideseq pipeline

"""

import os
import sys
import yaml
import argparse
import traceback
import subprocess
# Set up logger
import log
logger = log.createCustomLogger('root',"guideseq_V2")

from alignReads import alignReads
from filterBackgroundSites import filterBackgroundSites,filterBlackList,filterControl
from visualization import visualizeOfftargets
import identifyOfftargetSites
import validation
from tabulate import tabulate
from demultiplex import demultiplex,reformat_fastq

DEFAULT_YAML = os.path.dirname(os.path.realpath(__file__)) + "/default.yaml"

def get_parameters(manifest_data):

	# init
	default_refseqName = os.path.dirname(os.path.realpath(__file__)) + "/refseq_gene_name.py"
	with open(DEFAULT_YAML, 'r') as f:
		default = yaml.safe_load(f)
	with open(manifest_data, 'r') as f:
		return_dict = yaml.safe_load(f) # this is user input YAML
	default['analysis_folder'] = os.getcwd()
	default['refseq_names'] = default_refseqName
	default['Manhattan.R'] = os.path.dirname(os.path.realpath(__file__)) + "/Manhattan.R"
	validation.validateManifest(return_dict)
	
	# assign default if user YAML missing some parameters
	return_dict = assign_default(return_dict,default)
	return_dict['control_coord_chr'],tmp = return_dict['control_primer_coord'].split(":")
	return_dict['control_coord_start'],return_dict['control_coord_end'] = [int(x) for x in tmp.split("-")]

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

def get_I1_I2(R1,R2):
	"""Infer I1 I2 file name
	
	Assumption, the only difference between R1 R2 should be just "1" or "2"
	
	"""
	R1_list = list(R1)
	R2_list = list(R2)
	I1_list = list(R1)
	I2_list = list(R2)
	# can't directly mutation a string, have to be a list
	for i in range(len(R1_list)):
		if R1_list[i]!=R2_list[i]:
			if R1_list[i-1]=="r":
				I1_list[i-1]="i"
				I2_list[i-1]="i"
			else:
				I1_list[i-1]="I"
				I2_list[i-1]="I"
			return "".join(I1_list),"".join(I2_list)
		


class GuideSeq:

	def __init__(self):
		pass

	def parseManifest(self, manifest_path,sample='all',user_dict=None):
		logger.info('Loading manifest...')

		try:
			manifest_data = get_parameters(manifest_path)
			self.samples = {}
			self.parameters = manifest_data
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
		is_file_flag = os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['forward'])
		is_file_flag = os.path.isfile(is_file_flag)
		is_2_Undetermined_samples = False
		try:
			if os.path.isfile(os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['controlforward'])):
				is_2_Undetermined_samples = True
		except:
			is_2_Undetermined_samples = False
		try:
			
			if is_file_flag:
				if is_2_Undetermined_samples:
					self.demultiplex2()
				else:
					logger.info('Demultiplexing Undetermined files...')
					demultiplex(os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['forward']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['reverse']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['index1']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['index2']),
						barcode_dict,
						out_dir=self.parameters['demultiplex']['out_dir'],
						mismatch=self.parameters['demultiplex']['mismatch'],
						min_reads=self.parameters['demultiplex']['min_reads'],subsample_reads=self.parameters['demultiplex']['subsample_reads'])
					logger.info('Successfully demultiplexed reads.')
			else:
				logger.info("Undetermined file not found: %s"%(self.parameters['demultiplex']['forward']))
				logger.info("Undetermined file not found: %s"%(self.parameters['demultiplex']['reverse']))
				logger.info("Assuming individual fastq files are specified in barcode1, barcode2, controlbarcode1, and controlbarcode2")
				self.add_UMI_to_fastq_name()
		except Exception as e:
			logger.error('Error demultiplexing reads.')
			logger.error(traceback.format_exc())
			quit()

	def demultiplex2(self):
		"""demultiplex two Undetermined files, one for treatment, one for control
		
		
		"""
		barcode_dict = {}
		control_barcode_dict = {}
		for sample in self.samples:
			barcode1 = self.samples[sample]['barcode1']
			barcode2 = self.samples[sample]['barcode2']
			barcode_dict[sample] = [barcode1,barcode2]
			
			barcode1 = self.samples[sample]['controlbarcode1']
			barcode2 = self.samples[sample]['controlbarcode2']
			control_barcode_dict["Control_"+sample] = [barcode1,barcode2]
		try:
			logger.info('Demultiplexing 2 Undetermined files...')

			demultiplex(os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['forward']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['reverse']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['index1']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['index2']),
						barcode_dict,
						out_dir=self.parameters['demultiplex']['out_dir'],
						mismatch=self.parameters['demultiplex']['mismatch'],
						min_reads=self.parameters['demultiplex']['min_reads'],subsample_reads=self.parameters['demultiplex']['subsample_reads'])
			demultiplex(os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['controlforward']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['controlreverse']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['controlindex1']),
						os.path.join(self.parameters['demultiplex']['input_dir'],self.parameters['demultiplex']['controlindex2']),
						control_barcode_dict,
						out_dir=self.parameters['demultiplex']['out_dir'],
						mismatch=self.parameters['demultiplex']['mismatch'],
						min_reads=self.parameters['demultiplex']['min_reads'],subsample_reads=self.parameters['demultiplex']['subsample_reads'])
			logger.info('Successfully demultiplexed reads.')

		except Exception as e:
			logger.error('Error demultiplexing reads.')
			logger.error(traceback.format_exc())
			quit()



	def add_UMI_to_fastq_name(self):
		"""
		
		This step is done in demultiplex function when user provide Undertermined file.
		
		If user's barcode1/2 and controlbarcode1/2 are files, then we suppose that user
		has done demultiplexing and just provided individual fastq files
		"""
		barcode_dict = {}

		try:
			logger.info('Reformating fastq files...The output files will be in the demultiplex folder')
			for sample in self.samples:
				R1 = self.samples[sample]['barcode1']
				R2 = self.samples[sample]['barcode2']
				I1,I2 = get_I1_I2(R1,R2)
				reformat_fastq(R1, R2, I1, I2, sample, out_dir=self.parameters['demultiplex']['out_dir'],subsample_reads=self.parameters['demultiplex']['subsample_reads'])
				I1,I2 = get_I1_I2(R1,R2)
				R1 = self.samples[sample]['controlbarcode1']
				R2 = self.samples[sample]['controlbarcode2']
				reformat_fastq(R1, R2, I1, I2, "Control_"+sample, out_dir=self.parameters['demultiplex']['out_dir'],subsample_reads=self.parameters['demultiplex']['subsample_reads'])
			logger.info('Successfully reformat input reads.')

		except Exception as e:
			logger.error('Error reformating reads.')
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
					
					sample_alignment_path = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.dedup.sam')
					# (HG19_path, read1, read2, outfile,njobs=6,umi_tools="umi_tools",samtools="samtools",bwa="bwa"):
					alignReads(self.parameters['reference_genome'],
							   self.samples[sample]['read1'],
							   self.samples[sample]['read2'],
							   sample_alignment_path,
							   njobs=self.parameters['njobs'],
							   bwa=self.parameters['bwa'],
							   samtools=self.parameters['samtools'],
							   umi_tools=self.parameters['umi_tools']
							   )
					self.samples[sample]['aligned'] = sample_alignment_path
					# sample = "control_"+sample
					sample_alignment_path = os.path.join(self.parameters['analysis_folder'], 'aligned', "Control_"+sample + '.dedup.sam')
					alignReads(self.parameters['reference_genome'],
							   self.samples[sample]['controlread1'],
							   self.samples[sample]['controlread2'],
							   sample_alignment_path,
							   njobs=self.parameters['njobs'],
							   bwa=self.parameters['bwa'],
							   samtools=self.parameters['samtools'],
							   umi_tools=self.parameters['umi_tools']
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
						self.samples[sample]['aligned'] = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.dedup.sam')
						self.samples[sample]['controlaligned'] = os.path.join(self.parameters['analysis_folder'], 'aligned', "Control_"+sample + '.dedup.sam')
					# Prepare sample annotations
					sample_data = self.samples[sample]
					annotations = {}
					annotations['Description'] = sample_data['description']
					annotations['Targetsite'] = sample
					# print ("Using control primer",sample_data['control_primer'])

					annotations['Sequence'] = sample_data['target']
					# print (annotations)


					self.identified[sample] = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.txt')

					identifyOfftargetSites.analyze(self.samples[sample]['aligned'], self.parameters['reference_genome'], self.identified[sample], annotations,
												   self.parameters['window_size'], self.parameters['max_score'], self.parameters['control_primer'],self.parameters)

					sample = "Control_"+sample

					self.identified[sample] = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.txt')

					identifyOfftargetSites.analyze(self.samples[sample.replace("Control_","")]['controlaligned'], self.parameters['reference_genome'], self.identified[sample], annotations,
												   self.parameters['window_size'], self.parameters['max_score'], self.parameters['control_primer'],self.parameters)
				except:
					logger.error(f'Failed for sample {sample}')   
					logger.error(traceback.format_exc())

			logger.info('Finished identifying offtarget sites.')

		except Exception as e:
			logger.error('Error identifying offtarget sites.')
			logger.error(traceback.format_exc())
			quit()

	def filterBackgroundSites(self):
		logger.info('Filtering background sites')

		# self.filtered = {}
		for sample in self.samples:
			try:
				out = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
				filterBlackList(self.parameters['bedtools'], self.identified[sample], self.parameters['blacklist'], out)
				self.identified[sample] = out
				sample = "Control_"+sample
				out = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
				filterBlackList(self.parameters['bedtools'], self.identified[sample], self.parameters['blacklist'], out)
				self.identified[sample] = out
				
				# combine treatment and control
				sample = sample.replace("Control_","")
				out = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.with_control_counts.txt')
				filterControl(self.parameters['bedtools'], self.identified[sample], self.identified["Control_"+sample], out)
				logger.info('Finished background filtering for {0} sample'.format(sample))

			except Exception as e:
				logger.error('Error filtering background sites: %s'%(sample))
				logger.error(traceback.format_exc())

	def visualize(self):
		logger.info('Visualizing off-target sites')


		for sample in self.samples: ## 3/6/2020 Yichao solved: visualization stopped when one sample failed
			try:
				infile = self.identified[sample]
				outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '_offtargets')
				try:
					self.parameters['PAM']
					visualizeOfftargets(infile, outfile, title=sample,PAM=self.parameters['PAM'],genome=self.parameters['genome'],refseq_names=self.parameters['refseq_names'])
				except:
					visualizeOfftargets(infile, outfile, title=sample,PAM="NGG",genome=self.parameters['genome'],refseq_names=self.parameters['refseq_names'])
				# Manhattan plot
				outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '.Manhattan.pdf')
				# command = f"{self.parameters['Rscript']} {self.parameters['Manhattan.R']} {self.parameters['on_target_reference_sequence']} {infile.replace('.txt','.rm_no_match.annot.tsv')} {outfile}"
				command = f"{self.parameters['Rscript']} {self.parameters['Manhattan.R']} {self.samples[sample]['target']} {infile.replace('.txt','.rm_no_match.annot.tsv')} {outfile}"
				logger.info(command)
				subprocess.call(command,shell=True)
			except Exception as e:
				logger.error('Error visualizing off-target sites: %s'%(sample))
				logger.error(traceback.format_exc())
			try:
				sample = "Control_"+sample
				infile = self.identified[sample]
				outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '_offtargets')
				try:
					self.parameters['PAM']
					visualizeOfftargets(infile, outfile, title=sample,PAM=self.parameters['PAM'],genome=self.parameters['genome'],refseq_names=self.parameters['refseq_names'])
				except:
					visualizeOfftargets(infile, outfile, title=sample,PAM="NGG",genome=self.parameters['genome'],refseq_names=self.parameters['refseq_names'])
				# Manhattan plot
				outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '.Manhattan.pdf')
				# command = f"{self.parameters['Rscript']} {self.parameters['Manhattan.R']} {self.parameters['on_target_reference_sequence']} {infile.replace('.txt','.rm_no_match.annot.tsv')} {outfile}"
				command = f"{self.parameters['Rscript']} {self.parameters['Manhattan.R']} {self.samples[sample]['target']} {infile.replace('.txt','.rm_no_match.annot.tsv')} {outfile}"
				logger.info(command)
				subprocess.call(command,shell=True)
			except Exception as e:
				logger.error('Error visualizing off-target sites: %s'%(sample))
				logger.error(traceback.format_exc())
		logger.info('Finished visualizing off-target sites')
	def parallel(self, manifest_path, lsf,step,overwrite):
		logger.info('Submitting parallel jobs for GuideSeqV2')
		if not os.path.exists("HPC_parallel_log"):
			os.makedirs("HPC_parallel_log")
		current_script = __file__
		count = 1
		try:
			for sample in self.samples:
				cmd = f'python {current_script} main --manifest {manifest_path} --sample {sample} --step {step} --overwrite "{overwrite}"'
				# logger.info(cmd)
				# continue
				# subprocess.call(lsf.replace("{Sample_Name}",sample).split() + [f"-J {sample[:20]}"] + [cmd])
				cmd=lsf.replace("{Sample_Name}",sample) + f" -J {sample[:10]} " + cmd
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
	parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -n 12 -R "span[hosts=1] rusage[mem=12000]" -P GUIDEseqV2 -q standard -o HPC_parallel_log/GUIDEseqV2_{Sample_Name}_%J.log')
	parallel_parser.add_argument('--step', help='Specify which steps of pipepline to run (demultiplex, align, identify,visualize)', default='demultiplex+align+identify+visualize')
	parallel_parser.add_argument('--overwrite', help='overwrite specifications in the yaml file', default=None,type=yaml.load)

	main_parser = subparsers.add_parser('main', help='Run a single step or a series of steps for one or all samples')
	main_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	main_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')
	main_parser.add_argument('--step', help='Specify steps, demultiplex, align, identify,visualize, order does not matter', default='demultiplex+align+identify+visualize')
	main_parser.add_argument('--overwrite', help='overwrite specifications in the yaml file', default=None,type=yaml.load)
	
	init_parser = subparsers.add_parser('init', help='initialize a default yaml in the current folder')

	return parser.parse_args()


def main():
	args = parse_args()
	

	if args.command == "init":
		subprocess.call(f"cp {DEFAULT_YAML} Guideseq.yaml",shell=True)
		print("Guideseq.yaml has been initialize. \nPlease add sample info and then run GUIDE-seq analysis using the 'main' or the 'parallel' subcommand.")
	elif args.command == "main":
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
			g.filterBackgroundSites()
		if "visualize" in steps:
			try:
				g.identified
			except:
				g.identified = {}
				for sample in g.samples:
					if 'identified' in g.samples[sample]:
						g.identified[sample] = g.samples[sample]['identified']
						g.identified["Control_"+sample] = g.samples[sample]['controlidentified']
					else:
						g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
						sample = "Control_"+sample
						g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
			g.visualize()

	elif args.command == 'parallel':
		c = GuideSeq()
		c.parseManifest(args.manifest,args.sample,args.overwrite)
		# exit()
		c.parallel(args.manifest, args.lsf, args.step,args.overwrite)

 
if __name__ == '__main__':
	main()
