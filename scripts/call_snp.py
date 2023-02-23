"""
This script is to find the biased heterozygous sites between three mapping approaches
Output:


function:

snp_call_crude(reference_file, index_file):
snp_call_precise(reference_file, index_file,threshold):
nestedDict_to_file(nested_dict,index_file,output_name):

find_genes():
"""
from __future__ import division
import sys
import numpy
from Bio import SeqIO
import gzip
import pandas as pd

#reference_file = sys.argv[1]
#index_file = sys.argv[2]

def read_ref(reference_file):
	with gzip.open(reference_file, "rt") as rf:
		chrom_dict = {rec.id : rec.seq for rec in SeqIO.parse(rf, "fasta")}
	print("reference copied!")
	return chrom_dict

def snp_call_crude(chrom_dict, index_file):
	"""this function is meant to find the sites coordinations with considerable 
	essentially still the differences finding but this time with quality control (2/8)
	so still calling the difference -> check if the one of the samples are good quality heterozygousity """

	result_report = {}

	#set index name for comparison name between mapping approaches 
	#compare_index = ['bwa-aln_vs_vg-all','bwa-aln_vs_vg-mammoth','vg-all_vs_vg-mammoth']

	#get mammoth names
	mammoth_namelist = []
	with open(index_file,'r') as id:
		mammoth_name = ""
		for line in id:
			line = line.strip()
			mammoth_name = line.split("\t")[2]
			if (mammoth_name in mammoth_namelist) == False:
				mammoth_namelist.append(mammoth_name)

	for i in range(len(mammoth_namelist)):
		#mammoths were numbered as 0-3; 0,1,2 in the keys are corresponds to comparision in the compare_index
		result_report[i] = {0:[],1:[],2:[]} #0,1,2 corresponds to compare_index 

	#set names of chromosomes
	chr_list = []
	for i in range(1,28):
		chr_list.append('chr'+str(i))
	chr_list.append('chrX')

	for chromosome in chr_list:
		mpfile = "merged.{chrname}.mpile.gz".format(chrname = chromosome)
		with gzip.open(mpfile, "r") as f1:
			for line in f1:
				reference = ""
				line = line.decode()
				line = line.strip()
				splited = line.split("\t")
				reference = chrom_dict[chromosome][int(splited[1])-1]

				for i in range(2,2+len(mammoth_namelist)): #iterate through mammoths
					bwa_hetero_rate = 0
					vg_all_hetero_rate = 0
					vg_mm_hetero_rate = 0

					#ensure all samples are non-empty
					if (splited[i] != '*') and (splited[i+4] != '*') and (splited[i+8] != '*'):

						#bwa-aln vs vg_all (0 in compare_index)
						#first ensure that two mapping approaches have different results (different nucleotide identities)
						if set(splited[i]) != set(splited[i+4]):
							bwa_hetero_rate = round(splited[i].count(reference)/ len(splited[i]),4)
							vg_all_hetero_rate = round(splited[i+4].count(reference) / len(splited[i+4]),4)

							#0.2-0.8 is the crude threhold to determine a good heterzygous site
							#either one mapping approach has a good heterzygous site then this site coordination will be recorded
							if (0.2 < bwa_hetero_rate < 0.8) or (0.2 < vg_all_hetero_rate < 0.8):
								result_report[i-2][0].append(str(splited[0]+'_'+str(splited[1])))

						#bwa-aln vs vg_mm (1 in compare_index)
						#ensure both samples are non-empty
						#two mapping approaches have different results 
						if set(splited[i]) != set(splited[i+8]):
							bwa_hetero_rate = round(splited[i].count(reference)/ len(splited[i]),4)
							vg_mm_hetero_rate = round(splited[i+8].count(reference) / len(splited[i+8]),4)

							if (0.2 < bwa_hetero_rate < 0.8) or (0.2 < vg_mm_hetero_rate < 0.8):
								result_report[i-2][1].append(str(splited[0]+'_'+str(splited[1])))

						#vg_all vs vg_mm (2 in compare_index)
						#two mapping approaches have different results 
						if set(splited[i+4]) != set(splited[i+8]):
							vg_all_hetero_rate = round(splited[i+4].count(reference) / len(splited[i+4]),4)
							vg_mm_hetero_rate = round(splited[i+8].count(reference) / len(splited[i+8]),4)

							if (0.2 < vg_all_hetero_rate < 0.8) or (0.2 < vg_mm_hetero_rate < 0.8):
								result_report[i-2][2].append(str(splited[0]+'_'+str(splited[1])))

		print('{chrom} SNP calling is done!'.format(chrom = chromosome))
	print('All SNP calling is done!')
	return result_report

def nestedDict_to_file(nested_dict,index_file,output_name):
	"""this function is intended to parse nested dictionary into a result files"""
	compare_index = ['bwa-aln_vs_vg-all','bwa-aln_vs_vg-mammoth','vg-all_vs_vg-mammoth']

	mammoth_namelist = []
	with open(index_file,'r') as id:
		mammoth_name = ""
		for line in id:
			line = line.strip()
			mammoth_name = line.split("\t")[2]
			if (mammoth_name in mammoth_namelist) == False:
				mammoth_namelist.append(mammoth_name)

	with open(output_name,'w+') as w1:
		for mammoth_num in nested_dict: #for each mammoths 
			sample_name = ""
			for comparison_num, coordinates_list in nested_dict[mammoth_num].items():
				sample_name = mammoth_namelist[mammoth_num] + '_' + compare_index[comparison_num]
				w1.write(sample_name +'\t' + '\t'.join(coordinates_list) + '\n')



chrom_dict = read_ref(reference_file)

result_nested_dict = snp_call_crude(chrom_dict,index_file)
print('Coordinates obtained!')

nestedDict_to_file(result_nested_dict,index_file,'called_SNPs')

#test if the nestedDict_to_file works as expected 
"""
sample_nested = {0:{0:['chr1_10000','chr1_10001'], 1:['chr2_20000'], 2:['chr3_30000']},
1:{0:['chr1_11111'], 1:['chr2_21111'], 2:['chr3_31111']}, 
2:{0:['chr1_12222'], 1:['chr2_22222'], 2:['chr3_32222']}, 
3:{0:['chr1_13333'], 1:['chr2_22222'], 2:['chr3_33333']}}
nestedDict_to_file(sample_nested,index_file,'Called Potential New Heterozygous SNPs')
"""




















