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
from umi import umitag, consolidate
from visualization import visualizeOfftargets
import identifyOfftargetSites
import validation
from tabulate import tabulate

DEFAULT_WINDOW_SIZE = 25
DEFAULT_MAX_SCORE = 7

CONSOLIDATE_MIN_QUAL = 15
CONSOLIDATE_MIN_FREQ = 0.9

def get_parameters(manifest_data):
	default_yaml = os.path.dirname(os.path.realpath(__file__)) + "/default.yaml"
	default_refseqName = os.path.dirname(os.path.realpath(__file__)) + "/refseq_gene_name.py"
	with open(default_yaml, 'r') as f:
		default = yaml.load(f)
	with open(manifest_data, 'r') as f:
		return_dict = yaml.load(f)
	default['analysis_folder'] = os.getcwd()
	default['refseq_names'] = default_refseqName
	default['Manhattan.R'] = os.path.dirname(os.path.realpath(__file__)) + "/Manhattan.R"
	validation.validateManifest(return_dict)
	return_dict['parameters'] = {}
	for p in default:
		if not p in return_dict:
			return_dict['parameters'][p] = default[p]
		elif p!= "samples":
			return_dict['parameters'][p] = return_dict[p]
	return_dict['parameters']['control_coord_chr'],tmp = return_dict['parameters']['control_primer_coord'].split(":")
	return_dict['parameters']['control_coord_start'],return_dict['parameters']['control_coord_end'] = [int(x) for x in tmp.split("-")]
	return return_dict


class GuideSeq:

	def __init__(self):
		pass

	def parseManifest(self, manifest_path,sample='all'):
		logger.info('Loading manifest...')

		try:
			manifest_data = get_parameters(manifest_path)
			self.samples = {}
			self.parameters = manifest_data['parameters']
			if sample != "all":
				self.samples[sample] = manifest_data['samples'][sample]
			else:
				self.samples = manifest_data['samples']
			logger.info("\n"+tabulate([(k,v) for k,v in self.parameters.items()])) 
			# print(tabulate([(k,v) for k,v in self.undemultiplexed.items()])) 
			logger.info("\n"+tabulate([(k,v) for k,v in self.samples[list(self.samples.keys())[0]].items()])) 

		except Exception as e:
			logger.error(
				'Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
			logger.error(traceback.format_exc())
			sys.exit()


	def umitag(self):
		logger.info('umitagging reads...')

		try:
			self.umitagged = {}
			for sample in self.samples:
				try:
					self.umitagged[sample] = {}
					self.umitagged[sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'umitagged', sample + '.r1.umitagged.fastq')
					self.umitagged[sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'umitagged', sample + '.r2.umitagged.fastq')

					umitag.umitag(self.samples[sample]['read1'],
								  self.samples[sample]['read2'],
								  self.samples[sample]['index1'],
								  self.samples[sample]['index2'],
								  self.umitagged[sample]['read1'],
								  self.umitagged[sample]['read2'],
								  os.path.join(self.parameters['analysis_folder'], 'umitagged'))

					control_sample = "control_"+sample
					self.umitagged[control_sample] = {}
					self.umitagged[control_sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'umitagged', control_sample + '.r1.umitagged.fastq')
					self.umitagged[control_sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'umitagged', control_sample + '.r2.umitagged.fastq')

					umitag.umitag(self.samples[sample]['controlread1'],
								  self.samples[sample]['controlread2'],
								  self.samples[sample]['controlindex1'],
								  self.samples[sample]['controlindex2'],
								  self.umitagged[control_sample]['read1'],
								  self.umitagged[control_sample]['read2'],
								  os.path.join(self.parameters['analysis_folder'], 'umitagged'))
				except:
					logger.error(f'UMItag failed for sample {sample}')
					logger.error(traceback.format_exc())

			logger.info('Successfully umitagged reads.')
		except Exception as e:
			logger.error('Error umitagging')
			logger.error(traceback.format_exc())
			quit()

	def consolidate(self):
		logger.info('Consolidating reads...')

		try:
			self.consolidated = {}
			for sample in self.samples:
				try:
					self.consolidated[sample] = {}
					self.consolidated[sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'consolidated', sample + '.r1.consolidated.fastq')
					self.consolidated[sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'consolidated', sample + '.r2.consolidated.fastq')

					consolidate.consolidate(self.umitagged[sample]['read1'], self.consolidated[sample]['read1'], self.parameters['CONSOLIDATE_MIN_QUAL'], self.parameters['CONSOLIDATE_MIN_FREQ'])
					consolidate.consolidate(self.umitagged[sample]['read2'], self.consolidated[sample]['read2'], self.parameters['CONSOLIDATE_MIN_QUAL'], self.parameters['CONSOLIDATE_MIN_FREQ'])

					control_sample = "control_"+sample
					self.consolidated[control_sample] = {}
					self.consolidated[control_sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'consolidated', control_sample + '.r1.consolidated.fastq')
					self.consolidated[control_sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'consolidated', control_sample + '.r2.consolidated.fastq')

					consolidate.consolidate(self.umitagged[control_sample]['read1'], self.consolidated[control_sample]['read1'], self.parameters['CONSOLIDATE_MIN_QUAL'], self.parameters['CONSOLIDATE_MIN_FREQ'])
					consolidate.consolidate(self.umitagged[control_sample]['read2'], self.consolidated[control_sample]['read2'], self.parameters['CONSOLIDATE_MIN_QUAL'], self.parameters['CONSOLIDATE_MIN_FREQ'])

				except:
					logger.error(f'Consolidate failed for sample {sample}')
					logger.error(traceback.format_exc())

			logger.info('Successfully consolidated reads.')
		except Exception as e:
			logger.error('Error umitagging')
			logger.error(traceback.format_exc())
			quit()

	def alignReads(self):
		logger.info('Aligning reads...')

		try:
			self.aligned = {}
			for sample in self.samples:
				try:
					sample_alignment_path = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.sam')
					alignReads(self.parameters['bwa'],
							   self.parameters['reference_genome'],
							   self.consolidated[sample]['read1'],
							   self.consolidated[sample]['read2'],
							   sample_alignment_path,njobs=self.parameters['njobs'])
					self.aligned[sample] = sample_alignment_path
					sample = "control_"+sample
					sample_alignment_path = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.sam')
					alignReads(self.parameters['bwa'],
							   self.parameters['reference_genome'],
							   self.consolidated[sample]['read1'],
							   self.consolidated[sample]['read2'],
							   sample_alignment_path,njobs=self.parameters['njobs'])
					self.aligned[sample] = sample_alignment_path
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
					# Prepare sample annotations
					sample_data = self.samples[sample]
					annotations = {}
					annotations['Description'] = sample_data['description']
					annotations['Targetsite'] = sample
					# print ("Using control primer",sample_data['control_primer'])

					annotations['Sequence'] = sample_data['target']
					# print (annotations)

					samfile = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.sam')

					self.identified[sample] = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.txt')

					identifyOfftargetSites.analyze(samfile, self.parameters['reference_genome'], self.identified[sample], annotations,
												   self.parameters['window_size'], self.parameters['max_score'], sample_data['control_primer'],self.parameters)

					sample = "control_"+sample

					# print ("Using control primer",sample_data['control_primer'])


					samfile = os.path.join(self.parameters['analysis_folder'], 'aligned', 'control.sam')

					self.identified[sample] = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.txt')

					identifyOfftargetSites.analyze(samfile, self.parameters['reference_genome'], self.identified[sample], annotations,
												   self.parameters['window_size'], self.parameters['max_score'], sample_data['control_primer'],self.parameters)
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
				sample = "control_"+sample
				out = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
				filterBlackList(self.parameters['bedtools'], self.identified[sample], self.parameters['blacklist'], out)
				self.identified[sample] = out
				
				# combine treatment and control
				sample = sample.replace("control_","")
				out = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.with_control_counts.txt')
				filterControl(self.parameters['bedtools'], self.identified[sample], self.identified["control_"+sample], out)
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
				sample = "control_"+sample
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
	def parallel(self, manifest_path, lsf,step):
		logger.info('Submitting parallel jobs for GuideSeqV2')
		current_script = __file__
		count = 1
		try:
			for sample in self.samples:
				cmd = f'python {current_script} single --manifest {manifest_path} --sample {sample} --step {step}'
				logger.info(cmd)
				subprocess.call(lsf.split() + [f"-J {sample[:10]}"] + [cmd])
				count += 1
			logger.info('Finished job submission')

		except Exception as e:
			logger.error('Error submitting jobs.')
			logger.error(traceback.format_exc())


def parse_args():
	parser = argparse.ArgumentParser()

	subparsers = parser.add_subparsers(description='Individual Step Commands',
									   help='Use this to run individual steps of the pipeline',
									   dest='command')

	all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
	all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	parallel_parser = subparsers.add_parser('parallel', help='Run a single step or a series of steps in parallel, for each sample in the yaml file')
	parallel_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -R rusage[mem=60000] -P GUIDEV2 -q priority')
	parallel_parser.add_argument('--step', help='Specify which steps of pipepline to run (all, umitag, consolidate, align, identify,visualize)', default='umitag+consolidate+align+identify+visualize')

	single_parser = subparsers.add_parser('single', help='Run a single step or a series of steps (other than demultiplex) for a single sample or all samples')
	single_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	single_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')
	single_parser.add_argument('--step', help='Specify steps, umitag, consolidate, align, identify,visualize, order does not matter', default='umitag+consolidate+align+identify+visualize')

	# umitag_parser = subparsers.add_parser('umitag', help='UMI tag demultiplexed FASTQ files for consolidation')
	# all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	# all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	# consolidate_parser = subparsers.add_parser('consolidate', help='Consolidate UMI tagged FASTQs')
	# all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	# all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	# align_parser = subparsers.add_parser('align', help='Paired end read mapping to genome')
	# all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	# all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	# identify_parser = subparsers.add_parser('identify', help='Identify GUIDE-seq offtargets')
	# all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	# all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	# visualize_parser = subparsers.add_parser('visualize', help='Visualize off-target sites')
	# all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	# all_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')

	return parser.parse_args()


def main():
	args = parse_args()
	logger.info("User command: %s"%(" ".join(sys.argv)))

	if args.command == 'all':

		g = GuideSeq()
		g.parseManifest(args.manifest,args.sample)
		g.umitag()
		g.consolidate()
		g.alignReads()
		g.identifyOfftargetSites()
		g.filterBackgroundSites()
		g.visualize()


	elif args.command == "single":

		g = GuideSeq()
		g.parseManifest(args.manifest,args.sample)
		steps  = args.step.split("+")
		if "umitag" in steps:
			g.umitag()
		if "consolidate" in steps:
			g.consolidate()
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
					g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
					sample = "control_"+sample
					g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')


			g.visualize()

	elif args.command == 'parallel':
		c = GuideSeq()
		c.parseManifest(args.manifest)
		c.parallel(args.manifest, args.lsf, args.step)

 
if __name__ == '__main__':
	main()
