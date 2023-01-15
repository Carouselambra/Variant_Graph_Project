import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import ptitprince as pt


input_file = ""
#compare_dict = {}
def read_simple_file(input_file):
    """used to read file with sample name - single value pair"""
    compare_dict = {}
    with open(input_file, 'r') as f1:
        for line in f1:
            line = line.strip()
            splited = line.split('\t')
            compare_dict[splited[0]] = float(splited[1])

def read_to_df(filename):
    df = pd.read_csv(filename, index_col = 0)
    return df

def filter_dict(df):
    df = df.iloc[1:]
    df = df.iloc[ :-1]
    dict = df.to_dict()
    return dict

def simple_grouped_bar_plot_samples_together(compare_dict):
    mammoths = ['E467','L163','M6','oimyakon']
    bar1 = [v for k,v in compare_dict.items() if 'bwa-aln' in k]
    bar2 = [v for k,v in compare_dict.items() if 'variant-graph-all' in k]
    bar3 = [v for k,v in compare_dict.items() if 'variant-graph-mammoth-only' in k]

    barWidth = 0.25

    r1 = np.arange(len(bar1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]

    plt.bar(r1, bar1, color='#fac748', width=barWidth, edgecolor='white', label= 'bwa-aln')
    plt.bar(r2, bar2, color='#1d2f6f', width=barWidth, edgecolor='white', label='vg-all')
    plt.bar(r3, bar3, color='#8390fa', width=barWidth, edgecolor='white', label='vg-mammoth-only')

    plt.xlabel('mammoths', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bar1))], mammoths)

    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left",borderaxespad=0)

    plt.title("Whole Genome Nucleotide Identity Differences Across Three Mapping Strategies".format(type = type))

    plt.savefig('wg_vs_ref_nuc_id'.format(which = type),bbox_inches = "tight")
    plt.close()
    plt.clf()


def distribution_plot_by_mammoth(dict,data_name):
    mammoths = ['E467','L163','M6','oimyakon']
    for i in mammoths:
        #bar1 = [v for k,v in dict.items() if 'bwa-aln_'+i in k]
        bar1 = dict['bwa-aln_'+i].values()


        #bar2 = [v for k,v in dict.items() if 'variant-graph-all_'+i in k]
        bar2 = dict['variant-graph-all_'+i].values()

        #bar3 = [v for k,v in dict.items() if 'variant-graph-mammoth-only_'+i in k]
        bar3 = dict['variant-graph-mammoth-only_'+i].values()

        barWidth = 0.02

        r2 = dict['variant-graph-all_'+i].keys()

        r1 = [x - barWidth for x in r2]
        r3 = [x + barWidth for x in r2]

        plt.bar(r1, bar1, color='#fac748', width=barWidth, edgecolor='white', label= 'bwa-aln')
        plt.bar(r2, bar2, color='#1d2f6f', width=barWidth, edgecolor='white', label='vg-all')
        plt.bar(r3, bar3, color='#8390fa', width=barWidth, edgecolor='white', label='vg-mammoth-only')

        #plt.xticks(dict['variant-graph-all_'+i].keys())

        plt.legend()
        title_font = {'fontname':'Helvetica'}
        plt.title(i+' similarity rate distribution on non empty sites of three mapping approaches',title_font)

        plt.savefig('plot_zoomin_nofilter_'+i + '_' + data_name,bbox_inches='tight')
        plt.show()
        plt.cla()






    #ax.set_xticks([], [])


    plt.show()

def distribution_plot_by_mammoth_depth_series(filter_list):
    """this function is used to generate plots of heterzygousity distribution 
    by mammoths and series are filter changes"""

    #here have to read several depth files (different plot structure!)


data_name = "compare_rf_separate_hetero_distr_10bins_2nd_round"
distr_rate_df = read_to_df(data_name)
distribution_plot_by_mammoth(filter_dict(distr_rate_df),data_name.partition('_2nd')[0])


#distr_df = read_to_df('compare_ref_distribution_10bin')
#print(filter_dict(distr_df))
#distribution_plot_by_mammoth(filter_dict(distr_df))
"""
for i in range(4,24,2):
    data_name = 'compare_rf_sep_hetero_distr_10bins_deeper_than_{depth}'.format(depth = str(i))
    distr_df = read_to_df(data_name)
    distribution_plot_by_mammoth(filter_dict(distr_df),data_name)

    #print(filter_dict(distr_rate_df))

"""
#nested_bar_plot(filter_df(distr_df))
#print(distr_df)
