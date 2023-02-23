import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""This script is used to visualize the first results (coverage, nuecleotide difference rate, empty rates)

Input:
    all coverage file
    all nucleotide difference file
    all empty rate file

Out:
    plots
"""

#all_coverage_df = pd.read_csv('all_coverage',sep = '\t', index_col = 0)
#print(all_coverage_df)

#all_nuc_diff_df = pd.read_csv('all_nuc_diff',sep = '\t', index_col = 0)
#all_emptyrate_df = pd.read_csv('all_emptyrate',sep = '\t', index_col = 0)

# for sets with 4 mammoths * 3 mapping = 12 columns

def grouped_bar_plot_by_chromosomes(df,type):
    #generate 4 plots
    mammoths = ['E467','L163','M6','oimyakon']
    for i in range(0,len(mammoths)):
        bar1 = df.iloc[:,i] #bwa alignment
        bar2 = df.iloc[:,i+4] #vg-all
        bar3 = df.iloc[:,i+4+4] #vg-mammoth-only

        barWidth = 0.25

        r1 = np.arange(len(bar1))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]

        plt.bar(r1, bar1, color='#ff6f69', width=barWidth, edgecolor='white', label= 'bwa-aln')
        plt.bar(r2, bar2, color='#ffcc5c', width=barWidth, edgecolor='white', label='vg-all')
        plt.bar(r3, bar3, color='#88d8b0', width=barWidth, edgecolor='white', label='vg-mammoth-only')

        plt.xlabel('chromosomes', fontweight='bold')
        plt.xticks([r + barWidth for r in range(len(bar1))], list(df.index))

        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left",borderaxespad=0)

        plt.savefig('{type}_{mammoth}_results'.format(mammoth = mammoths[i], type = type),bbox_inches = "tight")
        plt.close()
        plt.clf()

def read_wg_file(filepath):
    wg_dict = {}
    with open(filepath, 'r') as f1:
        for line in f1:
            line = line.strip()
            splited = line.split('\t')
            wg_dict[splited[0]] = float(splited[1])
    return wg_dict

def grouped_bar_plot_wg(wg_dict,type):
    #generate one plots with 3 var
    mammoths = ['E467','L163','M6','oimyakon']
    bar_1 = [v for k,v in wg_dict.items() if 'bwa-aln' in k]
    bar_2 = [v for k,v in wg_dict.items() if 'variant-graph-all' in k]
    bar_3 = [v for k,v in wg_dict.items() if 'variant-graph-mammoth-only' in k]

    barWidth = 0.25

    r1 = np.arange(len(bar_1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    plt.bar(r1, bar_1, color='#ff6f69', width=barWidth, edgecolor='white', label= 'bwa-aln')
    plt.bar(r2, bar_2, color='#ffcc5c', width=barWidth, edgecolor='white', label='vg-all')
    plt.bar(r3, bar_3, color='#88d8b0', width=barWidth, edgecolor='white', label='vg-mammoth-only')

    plt.xlabel('mammoths', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bar_1))], mammoths)

    plt.ylabel('rate of reference without read mapping', fontweight = 'bold')

    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left",borderaxespad=0)

    plt.title("Whole Genome {type} of Three Mapping Strategies".format(type = type))

    plt.savefig('wg_{which}_results'.format(which = type),bbox_inches = "tight")
    plt.close()
    plt.clf()

def grouped_bar_plot_wg_diff(wg_nuc_diff,type):

    mammoths = ['E467','L163','M6','oimyakon']
    bar_1 = [v for k,v in wg_nuc_diff.items() if 'bwa-aln_vs_vg-all' in k]
    bar_2 = [v for k,v in wg_nuc_diff.items() if 'bwa-aln_vs_vg-mammoth' in k]
    bar_3 = [v for k,v in wg_nuc_diff.items() if 'vg-all_vs_vg-mammoth' in k]

    barWidth = 0.25

    r1 = np.arange(len(bar_1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    plt.bar(r1, bar_1, color='#ba274a', width=barWidth, edgecolor='white', label= 'bwa-aln_vs_vg-all')
    plt.bar(r2, bar_2, color='#841c26', width=barWidth, edgecolor='white', label='bwa-aln_vs_vg-mammoth')
    plt.bar(r3, bar_3, color='#2191fb', width=barWidth, edgecolor='white', label='vg-all_vs_vg-mammoth')

    plt.xlabel('mammoths', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bar_1))], mammoths)

    plt.legend()#bbox_to_anchor=(1.04, 1), loc="upper left",borderaxespad=0)
    plt.title('Whole Genome Nucleotide Difference {emptyornot} Between Three Mapping Approaches'.format(emptyornot = type))

    plt.savefig('wg_Nuc_diff_results_{which}'.format(which = type),bbox_inches = "tight")
    plt.close()
    plt.clf()



#wg_coverage = read_wg_file("whole_genome_coverage")
#grouped_bar_plot_wg(wg_coverage,"Average Depth")

#wg_empty = read_wg_file("whole_genome_emptyrate")
#grouped_bar_plot_wg(wg_empty,"rate of reference without read mapping")

#wg_nuc_diff = read_wg_file("whole_genome_nuc_diff")
#grouped_bar_plot_wg_diff(wg_nuc_diff,"include empty")

wg_nuc_diff_no_0 = read_wg_file("whole_genome_nuc_diff_exclude_empty")

#grouped_bar_plot_wg_diff(wg_nuc_diff_no_0,'At Non-empty Sites')
#grouped_bar_plot_by_chromosomes(all_coverage_df,"mean_coverage")
#grouped_bar_plot_by_chromosomes(all_emptyrate_df,'empty_rate')
