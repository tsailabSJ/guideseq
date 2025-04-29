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
from demultiplex import demultiplex

DEFAULT_YAML = os.path.dirname(os.path.realpath(__file__)) + "/default.yaml"

def get_parameters(manifest_data):

	# init
	default_refseqName = os.path.dirname(os.path.realpath(__file__)) + "/refseq_gene_name.py"
	with open(DEFAULT_YAML, 'r') as f:
		default = yaml.load(f)
	with open(manifest_data, 'r') as f:
		return_dict = yaml.load(f) # this is user input YAML
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
		try:
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
				command = f"{self.parameters['Rscript']} {self.parameters['Manhattan.R']} {self.parameters['on_target_reference_sequence']} {infile.replace('.txt','.rm_no_match.annot.tsv')} {outfile}"
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
				command = f"{self.parameters['Rscript']} {self.parameters['Manhattan.R']} {self.parameters['on_target_reference_sequence']} {infile.replace('.txt','.rm_no_match.annot.tsv')} {outfile}"
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
				logger.info(cmd)
				# continue
				subprocess.call(lsf.replace("{Sample_Name}",sample).split() + [f"-J {sample[:10]}"] + [cmd])
				count += 1
			logger.info('Finished job submission')

		except Exception as e:
			logger.error('Error submitting jobs.')
			logger.error(traceback.format_exc())


# yaml file
# dedup sam file
# sample name
# target sequence
# output name
annotations = {}
annotations['Description'] = "Description"
annotations['Targetsite'] = sys.argv[3]
annotations['Sequence'] = sys.argv[4]

output = os.path.join(os.getcwd(), sys.argv[5] + '.identified.tsv')

g = GuideSeq()
g.parseManifest(sys.argv[1])
identifyOfftargetSites.analyze(sys.argv[2], g.parameters['reference_genome'], output, annotations,g.parameters['window_size'], g.parameters['max_score'], g.parameters['control_primer'],g.parameters)


