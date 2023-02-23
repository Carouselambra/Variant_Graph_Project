"""
This script is intended to compare mapping results of each sample with African elephant reference genome

Input:
	index file
	reference file
	mpfiles (specified in separate functions)

Output:
	a file containing difference rates between reference and each samples (total counts)
	note: this time instead of boolean, we count the identity difference rates

Functions:
	
"""
from __future__ import division
import sys
import numpy
from Bio import SeqIO
import gzip
#import collections
import pandas as pd

reference_file = sys.argv[1]
index_file = sys.argv[2]

def read_ref(reference_file):
	with gzip.open(reference_file, "rt") as rf:
		chrom_dict = {rec.id : rec.seq for rec in SeqIO.parse(rf, "fasta")}
	print("reference copied!")
	return chrom_dict

def similarity_rate_simpler(reference_file, index_file):
	"""report average and distribution of similarity rate of all sample
	remember to change bin numbers"""
	avg_dict = {}
	distr_dict = {}
	key_name = []
	counter = 0

	#build dictionary for average similarity rate of all samples
	for i in range(2,14):
		avg_dict[i] = 0
		distr_dict[i] = {}
		for j in range(0,10):
			distr_dict[i][j/10] = 0

	#build dictionary for distribution of similarity rate

	#record sample names from index file
	with open(index_file,'r') as id:
		sample_name = ""
		for line in id:
			line = line.strip()
			sample_name = line.split("\t")[1] + "_" + line.split("\t")[2]
			key_name.append(sample_name)
	print("Index copied!")

	#store reference genome as a dictionary  (key - chromosome, value - sequence string)
	with gzip.open(reference_file, "rt") as rf:
		chrom_dict = {rec.id : rec.seq for rec in SeqIO.parse(rf, "fasta")}
	print("reference copied!")

	chr_list = []
	for i in range(1,28):
		chr_list.append('chr'+str(i))
	chr_list.append('chrX')

	#open each chromosome 
	for chr in chr_list: #iterate through each chromosome
		mpfile = "merged.{chrname}.mpile.gz".format(chrname = chr)
		#build dictionary for distribution of similarity rate
		with gzip.open(mpfile, "r") as f1:
			for line in f1:
				line = line.decode()
				line = line.strip()
				splited = line.split("\t")
				chr_pos = int(splited[1]) #locate the chromosome coordinate

				if chrom_dict[chr][chr_pos-1] != "N": #ensure reference location is not empty
					if any(sample == "*" for sample in splited[2:]) == False: #ensure all sample have coverage
						counter += 1
						for i in range(2,len(splited)):
							sim_count = splited[i].count(chrom_dict[chr][chr_pos-1]) # number of same nucleotide with reference
							sim_rate = round(sim_count / len(splited[i]),1) # rate of the same nucleotide / all nucleotide

							#add this similarity rate into the sum
							avg_dict[i] += sim_rate

							#count this similarity rate into its bin

							if sim_rate in distr_dict[i]:
								distr_dict[i][sim_rate] += 1


		print("{chr} is done!".format(chr = chr))

	#report average similarity rate in each sample
	rate_sum_list = []
	for sample in avg_dict:
		rate_sum_list.append(round(avg_dict[sample]/counter,5))
	rate_sum_dict = dict(zip(key_name,rate_sum_list))

	#report distribution in a nested dictionary: key - sample name, value - a freq dictionary

	report_distr_dict = {}

	for i in range(len(key_name)):
		report_distr_dict[key_name[i]] = distr_dict[i+2]


	return rate_sum_dict, report_distr_dict

#def similarity_rate_coverage_filtered():
	"""this function is to calculate similarity rate distribution for each samples with minimum coverage required"""

#def similarity_rate_avg_random_sample(chrom_dict,index_file):
	"""this function is to calculate average of similarity rate across shared read sites
	but with random sampling at every site!"""

def separate_hetero_distr(chrom_dict,index_file,bins):
	"""this function is to calculate the distribution of bias rate at heterozygous site for each samples 
	input:
		chrom_dict - a dictionary contain reference genome, key - chromosome, value - a string of nucleotide sequences
		index_file - the file path to the index file
		bins
			bin should be one of the following numbers: 10, 100, or 1000

	output:
		result_hetero_distr_dict
			a dictionary contains the occurences of all heterozygous rates (0-1) for all samples

	"""
	hetero_distr_dict = {}

	#determine the round digits by bins
	digits = len(str(bins))-1 

	for i in range(2,14):
		#each sample has a dictionary with similarity rates and 'counter' as keys, occurance of each rates as values
		hetero_distr_dict[i] = {}
		hetero_distr_dict[i]['counter'] = 0
		for j in range(0,bins+1):
			hetero_distr_dict[i][j/bins] = 0

	chr_list = []
	for i in range(1,28):
		chr_list.append('chr'+str(i))
	chr_list.append('chrX')

	key_name = []
	with open(index_file,'r') as id:
		sample_name = ""
		for line in id:
			line = line.strip()
			sample_name = line.split("\t")[1] + "_" + line.split("\t")[2]
			key_name.append(sample_name)

	for chr in chr_list:
		mpfile = "merged.{chrname}.mpile.gz".format(chrname = chr)
		with gzip.open(mpfile, "r") as f1:
			for line in f1:
				line = line.decode()
				line = line.strip()
				splited = line.split("\t")
				chr_pos = int(splited[1]) #locate the chromosome coordinates

				if chrom_dict[chr][chr_pos-1] != "N": #ensure reference location is not empty
					for i in range(2,len(splited)):
						sim_count = splited[i].count(chrom_dict[chr][chr_pos-1])
						hetero_rate = round(1 - (sim_count / len(splited[i])),digits)
						if hetero_rate in hetero_distr_dict[i]:
							hetero_distr_dict[i][hetero_rate] += 1
							hetero_distr_dict[i]['counter'] += 1

		print("{chromosome} is done!".format(chromosome = chr))

	result_hetero_distr_dict = {}
	for sample in hetero_distr_dict: #access each sample by numeric order (e.g. sample number 1 is bwa-aln_E467)
		result_hetero_distr_dict[key_name[sample-2]] = {}
		for rate, count in hetero_distr_dict[sample].items():
			if rate != 'counter':
				result_hetero_distr_dict[key_name[sample-2]][rate] = round(hetero_distr_dict[sample][rate] / hetero_distr_dict[sample]['counter'],4)

	return result_hetero_distr_dict


def write_to_file_sum(result_dict,type):
	with open("compare_ref_{which}".format(which = type),'w') as w:
		for k,v in result_dict.items():
			w.write(str(k) + '\t' + str(v) + '\n')

def write_to_file_distr(nested_dict, filename):
	df = pd.DataFrame(nested_dict)
	df.to_csv(filename)

def separate_hetero_distr_filter_by_depth(chrom_dict,index_file,bins,depth_low):
	"""this function is to calculate the distribution of bias rate at heterozygous site for each samples
	 with additional minimum depth filter 
	input:
		chrom_dict 
			a dictionary contain reference genome, key - chromosome, value - a string of nucleotide sequences
		index_file  
			the file path to the index file
		bins (int, 10/100/1000)
			bin should be one of the following numbers: 10, 100, or 1000
		depth_low (int)
			the minimum depth required to be considered as a good occurance count
	"""
	hetero_distr_dict = {}

	#determine the round digits by bins
	digits = len(str(bins))-1 

	for i in range(2,14):
	#each sample has a dictionary with similarity rates and 'counter' as keys, occurance of each rates as values
		hetero_distr_dict[i] = {}
		hetero_distr_dict[i]['counter'] = 0
		for j in range(0,bins+1):
			hetero_distr_dict[i][j/bins] = 0

	chr_list = []
	for i in range(1,28):
		chr_list.append('chr'+str(i))
	chr_list.append('chrX')

	key_name = []
	with open(index_file,'r') as id:
		sample_name = ""
		for line in id:
			line = line.strip()
			sample_name = line.split("\t")[1] + "_" + line.split("\t")[2]
			key_name.append(sample_name)

	for chr in chr_list:
		mpfile = "merged.{chrname}.mpile.gz".format(chrname = chr)
		with gzip.open(mpfile, "r") as f1:
			for line in f1:
				line = line.decode()
				line = line.strip()
				splited = line.split("\t")
				chr_pos = int(splited[1]) #locate the chromosome coordinate


				if chrom_dict[chr][chr_pos-1] != "N": #ensure reference location is not empty
					for i in range(2,len(splited)): 
						#ensure this sample at this location pass the occurance threshold 
						if len(splited[i]) >= depth_low:
							sim_count = splited[i].count(chrom_dict[chr][chr_pos-1])
							hetero_rate = round(1 - (sim_count / len(splited[i])),digits)
							if hetero_rate in hetero_distr_dict[i]:
								hetero_distr_dict[i][hetero_rate] += 1
								hetero_distr_dict[i]['counter'] += 1


		print("{chr} is done!".format(chr = chr))

	result_hetero_distr_dict = {}
	for sample in hetero_distr_dict:
		result_hetero_distr_dict[key_name[sample-2]] = {}
		for rate, count in hetero_distr_dict[sample].items():
			if rate != 'counter':
				result_hetero_distr_dict[key_name[sample-2]][rate] = round(hetero_distr_dict[sample][rate] / hetero_distr_dict[sample]['counter'],4)

	print('full result dict is obtained!')

	return result_hetero_distr_dict

def table_data_heterozygous_rate_count(file):
	df = pd.read_csv(file,index_col=0)
	for sample in df:
		sample_hetero_rate = df[sample].iloc[0] + df[sample].iloc[-1]
		print(1-sample_hetero_rate)

#table_data_heterozygous_rate_count('compare_rf_separate_hetero_distr_10bins_2nd_round')
table_data_heterozygous_rate_count('compare_rf_sep_hetero_distr_10bins_deeper_than_16')



"""
def write_to_file_all_rates(all_sim_rate_dict):
		df = pd.DataFrame.from_dict(all_sim_rate_dict, orient = 'columns')
		df.to_csv("compare_ref_all_sim")
"""

#sim_avg, sim_distr = similarity_rate_simpler(reference_file, index_file)

#write_to_file_sum(sim_avg, 'avg_sim_rate_simpler')
#write_to_file_distr(sim_distr)


#reference_dict = read_ref(reference_file)

#to obtain heterozygous rate distribution without depth filter
#NOTE: results here are ratio (occurance of heterozygous / occurance of total) !!!

#hetero_distr_10 = separate_hetero_distr(reference_dict,index_file,10)
#write_to_file_distr(hetero_distr_10,'compare_rf_separate_hetero_distr_10bins_2nd_round')

#hetero_distr_100 = separate_hetero_distr(reference_dict,index_file,100)
#write_to_file_distr(hetero_distr_100,'compare_rf_separate_hetero_distr_100bins_2nd_round')


#heterozygous rate distribution filtered by depth 
"""
depths_range = []
for i in range(14,24,2):
	depths_range.append(i)
for depth in depths_range:
	hetero_distr_filtered = separate_hetero_distr_filter_by_depth(reference_dict,index_file,10,depth)
	write_to_file_distr(hetero_distr_filtered,'compare_rf_sep_hetero_distr_10bins_deeper_than_{num}'.format(num = str(depth)))
	print('Minimum depth filter {num} is applied!'.format(num = depth))
"""
