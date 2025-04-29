import subprocess
import os


def filterBackgroundSites(bedtools_path, sample_path, control_path, outfile):
	output_folder = os.path.dirname(outfile)
	if not os.path.exists(output_folder):
		try:
			os.makedirs(output_folder)
		except:
			print ("Folder exist")

	bedtools_filter_command = '{0} intersect -a {1} -b {2}'.format(bedtools_path, sample_path, control_path)
	print (bedtools_filter_command)
	with open(outfile, 'w') as outfile:
		subprocess.call(bedtools_filter_command.split(), stdout=outfile)

def filterControl(bedtools_path, sample_path, control_path, outfile):
	output_folder = os.path.dirname(outfile)
	if not os.path.exists(output_folder):
		try:
			os.makedirs(output_folder)
		except:
			print ("Folder exist")

	bedtools_sort_command = '{0} sort -i {1} > {1}.sorted'.format(bedtools_path, sample_path)
	bedtools_sort_command1 = '{0} sort -i {1} > {1}.sorted'.format(bedtools_path, control_path)
	bedtools_closest_command = '{0} closest -a {1}.sorted -b {2}.sorted -d > {3};rm {1}.sorted;rm {2}.sorted'.format(bedtools_path, sample_path, control_path,outfile)

	subprocess.call(bedtools_sort_command,shell=True)
	subprocess.call(bedtools_sort_command1,shell=True)
	subprocess.call(bedtools_closest_command,shell=True)

def filterBlackList(bedtools_path, sample_path, bl_path, outfile):
	output_folder = os.path.dirname(outfile)
	if not os.path.exists(output_folder):
		try:
			os.makedirs(output_folder)
		except:
			print ("Folder exist")
	
	bedtools_filter_command = '{0} intersect -a {1} -b {2} -v -wa -header > {3}'.format(bedtools_path, sample_path, bl_path,outfile)
	# print (bedtools_filter_command)
	subprocess.call(bedtools_filter_command,shell=True)



