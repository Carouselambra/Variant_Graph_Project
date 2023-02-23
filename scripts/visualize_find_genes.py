from __future__ import division
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import ptitprince as pt

"""gene_counts_to_all
"""
def data_gene_counts_to_all_no_filter():
	mammoths = ['E467','L163','M6','oimyakon']
	compare_list = ['bwa-aln_vs_vg-all','bwa-aln_vs_vg-mammoth','vg-all_vs_vg-mammoth']

	total_genes = sum(1 for line in open('loxafr4_genes_added_mouse_orthologs.gff3'))
	gene_rate_dict = {}

	for m in mammoths:
		for c in compare_list:
			file_name = 'foundgenes_{mammoth}_{compare}'.format(mammoth = m, compare = c)
			sample_gene_count = sum(1 for line in open(file_name))
			sample_gene_rate = round(sample_gene_count/total_genes,4)
			gene_rate_dict['{mm}_{compare}'.format(mm = m, compare = c)] = sample_gene_rate
	return gene_rate_dict


def plot_gene_counts_to_all(gene_rate_dict):
	"""this is to plot the output of genes ratio"""
	mammoths = ['E467','L163','M6','oimyakon']
	bar_1 = [v for k,v in gene_rate_dict.items() if 'bwa-aln_vs_vg-all' in k]
	bar_2 = [v for k,v in gene_rate_dict.items() if 'bwa-aln_vs_vg-mammoth' in k]
	bar_3 = [v for k,v in gene_rate_dict.items() if 'vg-all_vs_vg-mammoth' in k]

	barWidth = 0.25

	r1 = np.arange(len(bar_1))
	r2 = [x + barWidth for x in r1]
	r3 = [x + barWidth for x in r2]

	plt.bar(r1, bar_1, color='#FFC759', width=barWidth, edgecolor='white', label= 'bwa-aln_vs_vg-all')
	plt.bar(r2, bar_2, color='#D2B10F', width=barWidth, edgecolor='white', label='bwa-aln_vs_vg-mammoth')
	plt.bar(r3, bar_3, color='#65AFFF', width=barWidth, edgecolor='white', label='vg-all_vs_vg-mammoth')

	plt.xlabel('mammoths', fontweight='bold')
	plt.xticks([r + barWidth for r in range(len(bar_1))], mammoths)

	plt.ylabel('percentage of genes associated')

	plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left",borderaxespad=0)
	plt.title('Whole Genome Genes Difference Between Three Mapping Approaches')

	plt.savefig('wg_genes_diff',bbox_inches = "tight")
	plt.close()
	plt.clf()



print(data_gene_counts_to_all_no_filter())
plot_gene_counts_to_all(data_gene_counts_to_all_no_filter())