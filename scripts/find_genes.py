"""This script is used to find genes with SNPs in eacn samples"""

from __future__ import division
import sys
import numpy
from Bio import SeqIO
import gzip
import pandas as pd

#called_SNPs_file = sys.argv[1]
#all_gene_file = sys.argv[2]

def extract_coordinates(SNPs_file):
	"""this function is to extract all SNPs found in call_snps.py """
	all_SNPs = {}
	with open(SNPs_file,'r') as f:
		for line in f:
			line = line.strip()
			sample_name = line.split('\t')[0]
			locs = line.split('\t')[1:]
			all_SNPs[sample_name] = locs

	return all_SNPs

def find_genes(all_SNPs,all_gene_file):
	"""this function is take extracted coordinates locate any genes in the region 
	use genes 
	output format:
	sample (mammonth-comparision) per file? 
	you'll worry about which file to look at later please.
	key - value pairs genes - SNPs
	"""
	with open(all_gene_file,'r') as f:
		for line in f:
			chrom = ""
			start_loc = 0
			end_loc = 0
			gene_info = ""

			line = line.strip()
			splited = line.split("\t")

			chrom = splited[0]
			start_loc = int(splited[3])
			end_loc = int(splited[4])
			if len(splited) >= 9:
				gene_info = str(splited[8]).split(";")
				gene_id = str(gene_info[0].partition(':')[2])
				gene_name = str(gene_info[1].partition('=')[2])
				gene_identifer = gene_id + '_' + gene_name
				gene_len = end_loc - start_loc + 1

			#iterate through _SNPs dict 
			
				for sample, coord_list in all_SNPs.items():
					#each individual sample has its own file 
					out_file = "full_foundgenes_{sample_n}".format(sample_n = sample)
					sample_gene_dict = {}
					for SNP in coord_list:
						if chrom == SNP.split('_')[0]: #match the searched gene with SNPs at chromosome
							if start_loc <= int(SNP.split('_')[1]) <= end_loc: #match the SNP is within the gene range
								if gene_identifer in sample_gene_dict:
									sample_gene_dict[gene_identifer].append(SNP)
								else:
									sample_gene_dict[gene_identifer] = []
									sample_gene_dict[gene_identifer].append(SNP)
								
								#print("got a hit on {this} for this {compare}".format(this = SNP, compare = sample))
								#entry = str(chrom) + '\t' + str(SNP.split('_')[1]) + '\t' + str(gene_info)

					with open(out_file,'a') as o_f:
						for gene_identifer, SNPs in sample_gene_dict.items():
							#index:gene name - gene id - gene length - SNPs counts - SNP/bp rate - SNPs coordinates 
							o_f.write(gene_name + '\t' + gene_id + '\t' + str(gene_len) +'\t' + str(len(SNPs)) + '\t' + str(round(len(SNPs)/gene_len,8)) + '\t' + ','.join(SNPs) + '\n')
									
def find_genes_report(foundgenes_files):
	for file in foundgenes_files:
		gene_num_count = sum(1 for line in open(file))

def check_genes_overview(foundgenes_file, min_SNP_per_bp):
	"""this funtion is used to report 
	how many genes are found in each sample has SNPs/bp counts above the set minimum
	"""
	col_names = ['gene_name','gene_id','gene_length','SNPs_counts','SNP_per_bp','SNPs']
	sample_df = pd.read_csv(foundgenes_file,names = col_names, sep = '\t')
	sample_df['SNPs_counts'] = pd.to_numeric(sample_df['SNPs_counts'])
	sample_df['gene_length'] = pd.to_numeric(sample_df['gene_length'])
	sample_df['SNP_per_bp'] = pd.to_numeric(sample_df['SNP_per_bp'])

	sorted_df = sample_df.sort_values(by = ['SNP_per_bp'], ascending= False)
	filtered_df = sorted_df.loc[sorted_df['SNP_per_bp'] >= min_SNP_per_bp]

	gene_num = len(filtered_df.index)
	print('{sample} has total {genes} genes with SNP per base pair more than {min}'.format(sample = foundgenes_file.partition('genes_')[2], genes = gene_num, min = min_SNP_per_bp))

def check_genes_top_highest_SNPperbp(foundgenes_file, topN):
	"""this function is used to find the genes where 
	- find the top N genes with most SNPs located in each samples
	"""
	col_names = ['gene_name','gene_id','gene_length','SNPs_counts','SNP_per_bp','SNPs']
	sample_df = pd.read_csv(foundgenes_file,names = col_names, sep = '\t')
	sample_df['SNPs_counts'] = pd.to_numeric(sample_df['SNPs_counts'])
	sample_df['gene_length'] = pd.to_numeric(sample_df['gene_length'])
	sample_df['SNP_per_bp'] = pd.to_numeric(sample_df['SNP_per_bp'])

	final_df = sample_df.sort_values(by = ['SNP_per_bp'], ascending= False)
	out_df = final_df.head(topN)

	o_f = 'Top{n}_{sample}'.format(n = topN, sample = foundgenes_file)
	out_df.to_csv(o_f,sep = '\t', index = False, header = False)

def top_genes_quick_report(topN_shared_files):
	"""quick report all gene names (simply add protein coding)"""	
	for file in topN_shared_files:
		name_list = []
		prtn_c_c = 0
		with open(file,'r') as f:
			for line in f:
				splited = line.split('\t')
				name = splited[0]
				if name == "protein_coding":
					prtn_c_c += 1
				else:
					name_list.append(name)
		out_entry = ', '.join(name_list) +', and {n} novel protein coding genes'.format(n = prtn_c_c)
		sample = file.partition('full_foundgenes_')[2]
		print(sample +": " +out_entry)




#all_SNPs = extract_coordinates(called_SNPs_file)
#find_genes(all_SNPs,all_gene_file)

mammoths = ['E467','L163','M6','oimyakon']
compare = ['bwa-aln_vs_vg-all','bwa-aln_vs_vg-mammoth', 'vg-all_vs_vg-mammoth']
foundgenes_files_list = []
topN_file_list = []
for m in mammoths:
	for c in compare:
		foundgenes_files_list.append('full_foundgenes_'+m+'_'+c)
		topN_file_list.append('Top20_full_foundgenes_'+m+'_'+c)

top_genes_quick_report(topN_file_list)
#for f in foundgenes_files_list:
	#check_genes_overview(f,0.1)
	#check_genes_top_highest_SNPperbp(f,20)

#check_genes_unique_in_vg_mm()


