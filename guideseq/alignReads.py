"""
alignReads
"""

import subprocess
import os
import logging
import sys

logger = logging.getLogger('root')
logger.propagate = False


def trim_reads(cutadapt=None, R1=None, R2=None, label=None,output_dir=None, njobs=None,error=0.1,overlap=20,dsODN_rev="GTTTAATTGAGTTGTCATATGTTAATAACGGTAT",dsODN_fwd="ATACCGTTATTAACATATGACAACTCAATTAAAC", **kwargs):
	"""Trimming PE reads and return trimmed file name"""
	logger.info(f'Start trimming for {label}')
	trimming_input = f"-G rev=^{dsODN_rev} -G fwd=^{dsODN_fwd} -e {error} --rename '{{r2.match_sequence}}_{{r2.adapter_name}}_{{id}}' --overlap {overlap}"
	trimmed_reads_PE_output = f"-o {output_dir}/{label}.R1.fastq --untrimmed-output {output_dir}/{label}.R1.untrimmed.fastq -p {output_dir}/{label}.R2.fastq --untrimmed-paired-output {output_dir}/{label}.R2.untrimmed.fastq"
	other_options = f"-j {njobs} --info-file {output_dir}/{label}.cutadapt_info_file.tsv"

	cutadapt_command = f"{cutadapt} {trimming_input} {trimmed_reads_PE_output} {other_options} {R1} {R2}"
	logger.info(cutadapt_command)
	# no_adapter
	subprocess.call(cutadapt_command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
	# exit()
	return f"{output_dir}/{label}.R1.fastq",f"{output_dir}/{label}.R2.fastq",f"{output_dir}/{label}.R1.untrimmed.fastq",f"{output_dir}/{label}.R2.untrimmed.fastq"

def alignReads(bwa=None, samtools=None, # Tools
			   reference_genome=None, R1=None, R2=None, label=None, output_dir=None, # Inputs
			   umi_tools=None, njobs=None,mapq_threshold=None, **kwargs): # other parameters
	"""Only for paired-end reads, return bam file name for next analysis

	Tools and inputs, e.g., reference_genome, R1, R2, can be absolute or relative path
	"""
	if not os.path.exists(output_dir):
		try:
			os.makedirs(output_dir)
			print ("creating folder",output_dir)
		except:
			print ("Folder exist")
	# Check if genome is already indexed by bwa
	index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']
	genome_indexed = True
	for extension in index_files_extensions:
		if not os.path.isfile(reference_genome + extension):
			genome_indexed = False
			break
	# If the genome is not already indexed, index it
	if not genome_indexed:
		logger.info('Genome index files not detected. Running BWA to generate indices.')
		bwa_index_command = f'{bwa} index {reference_genome}'
		logger.info(f'Running bwa command: {bwa_index_command}')
		subprocess.call(bwa_index_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
		logger.info('BWA genome index generated')

	# print (kwargs)
	R1,R2,noise_R1,noise_R2 = trim_reads(R1=R1, R2=R2, label=label, output_dir=output_dir,njobs=njobs, **kwargs)
	### alignment ###
	input_R1 = R1
	input_R2 = R2
	output_bam_file = f"{output_dir}/{label}.bam" # tmp
	output_sorted_bam_file = f"{output_dir}/{label}.st.bam"
	output_UMIdedup_bam_file = f"{output_dir}/{label}.dedup.bam"
	output_UMIdedup_sam_file = f"{output_dir}/{label}.dedup.sam"
	
	# BWA alignment GPU, 5/1 add MAPQ
	if bwa=="fq2bam":
		# bwa_alignment_command=f'pbrun fq2bam --ref {reference_genome} --in-fq {input_R1} {input_R2} --gpusort --num-gpus 1 --out-bam {output_sorted_bam_file}'
		# logger.info("GPU version BWA alignment")
		# logger.info(bwa_alignment_command)
		# subprocess.call(bwa_alignment_command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
		bwa_alignment_command=f'pbrun fq2bam --ref {reference_genome} --in-fq {input_R1} {input_R2} --gpusort --num-gpus 1 --out-bam {output_bam_file}'
		logger.info("GPU version BWA alignment")
		logger.info(bwa_alignment_command)
		subprocess.call(bwa_alignment_command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
		samtools_sort_command = f'{samtools} view -h -b -q {mapq_threshold} {output_bam_file} > {output_sorted_bam_file}'
		samtools_index_command = f"{samtools} index {output_sorted_bam_file}"
		logger.info('samtools applying MAPQ threshold')
		logger.info(samtools_sort_command)
		logger.info(samtools_index_command)
		subprocess.call(samtools_sort_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
		subprocess.call(samtools_index_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
		
	else:
		# BWA alignment dsODN reads
		opts = f"-t {njobs}"
		bwa_alignment_command = f'{bwa} mem {opts} {reference_genome} {input_R1} {input_R2} | {samtools} view -bS - > {output_bam_file}'
		logger.info("BWA alignment")
		logger.info(bwa_alignment_command)
		subprocess.call(bwa_alignment_command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
		
		# samtools sort
		samtools_sort_command = f'{samtools} view -h -b -q {mapq_threshold} {output_bam_file} | {samtools} sort -@ {njobs} -o {output_sorted_bam_file}'
		samtools_index_command = f"{samtools} index {output_sorted_bam_file}"
		logger.info('samtools sort bam file')
		logger.info(samtools_sort_command)
		logger.info(samtools_index_command)
		subprocess.call(samtools_sort_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
		subprocess.call(samtools_index_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)

	# UMI dedup
	UMItools_command = f"{umi_tools} dedup --stdin={output_sorted_bam_file} --paired --stdout={output_UMIdedup_bam_file} --method unique;{samtools} index {output_UMIdedup_bam_file};{samtools} view {output_UMIdedup_bam_file} > {output_UMIdedup_sam_file}"
	logger.info(UMItools_command)
	subprocess.call(UMItools_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	logger.info('Paired end mapping for {0} completed.'.format(label))


	### alignment ###
	input_R1 = noise_R1
	input_R2 = noise_R2
	label=f"{label}.mispriming"
	output_bam_file = f"{output_dir}/{label}.bam" # tmp
	output_sorted_bam_file = f"{output_dir}/{label}.st.bam"
	output_UMIdedup_bam_file = f"{output_dir}/{label}.dedup.bam"
	output_UMIdedup_sam_file = f"{output_dir}/{label}.dedup.sam"

	# BWA alignment GPU
	if bwa=="fq2bam":
		bwa_alignment_command=f'pbrun fq2bam --ref {reference_genome} --in-fq {input_R1} {input_R2} --gpusort --num-gpus 1 --out-bam {output_bam_file}'
		logger.info("GPU version BWA alignment")
		logger.info(bwa_alignment_command)
		subprocess.call(bwa_alignment_command,shell=True,stdout=sys.stdout,stderr=sys.stderr)
		samtools_sort_command = f'{samtools} view -h -b -q {mapq_threshold} {output_bam_file} > {output_sorted_bam_file}'
		samtools_index_command = f"{samtools} index {output_sorted_bam_file}"
		logger.info('samtools applying MAPQ threshold')
		logger.info(samtools_sort_command)
		logger.info(samtools_index_command)
		subprocess.call(samtools_sort_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
		subprocess.call(samtools_index_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	else:
		# BWA alignment dsODN reads
		opts = f"-t {njobs}"
		bwa_alignment_command = f'{bwa} mem {opts} {reference_genome} {input_R1} {input_R2} | {samtools} view -bS - > {output_bam_file}'
		logger.info("BWA alignment")
		logger.info(bwa_alignment_command)
		subprocess.call(bwa_alignment_command,shell=True,stdout=sys.stdout,stderr=sys.stderr)

		# samtools sort
		samtools_sort_command = f'{samtools} sort -@ {njobs} -o {output_sorted_bam_file} {output_bam_file}'
		samtools_index_command = f"{samtools} index {output_sorted_bam_file}"
		logger.info('samtools sort bam file')
		logger.info(samtools_sort_command)
		logger.info(samtools_index_command)
		subprocess.call(samtools_sort_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
		subprocess.call(samtools_index_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)

	# UMI dedup
	UMItools_command = f"{umi_tools} dedup --stdin={output_sorted_bam_file} --paired --stdout={output_UMIdedup_bam_file} --method unique;{samtools} index {output_UMIdedup_bam_file};{samtools} view {output_UMIdedup_bam_file} > {output_UMIdedup_sam_file}"
	logger.info(UMItools_command)
	subprocess.call(UMItools_command, shell=True,stdout=sys.stdout,stderr=sys.stderr)
	logger.info('Paired end mapping for {0} completed.'.format(label))


