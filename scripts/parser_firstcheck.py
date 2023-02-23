"""This script is used to check depth, empty mappings rates 
and nucleotide differences between mapping appraoches 
per chromosome in the mapping result files
Input:
    mpfile_file path
    index_file path

Output:
    coverage check file per chromosome (e.g. 'chr10_coverage_check.txt')
    nucleotide difference file per chromosome (e.g. 'chr10_nuc_diff.txt')
    empty mapping file per chromosome (e.g. 'chr10_emptyrate.txt')

Functions:
    read_file_to_df(mpfile_file, index_file)



"""
from __future__ import division
import sys
import pandas as pd
import numpy as np

#specify the input files 
mpfile_file = sys.argv[1]
index_file = sys.argv[2]

#extract name of chromosome 
id = mpfile_file.partition('merged.')[2].partition('.mpile')[0]


def read_file_to_df(mpfile_file, index_file):
    """This function is to used to convert a single mpfile to pandas dataframe
    input:
        mpfile_file: file path to mpfile
        index_file: file path to index file 
    output:
        all_vcf_df: a pandas dataframe with column names same with index file, contents from mpfile"""
    samples = ['chromosome','position']
    with open(index_file) as index:
        for line in index:
            line = line.strip()
            splited = line.split('\t')
            sample_name = '{mapping_method}_{mammoth}'.format(mapping_method = splited[1].strip(), mammoth = splited[2].strip())
            samples.append(sample_name)
    all_vcf_df = pd.read_csv(mpfile_file,compression = "gzip", sep = '\t', names = samples)
    print('VCF input done!')
    print(all_vcf_df)
    return all_vcf_df

def check_coverage_empty_rate(all_vcf_df):
    """This function is used to check average depth, and empty mapping rate of each chromosome
    input:
        all_vcf_df: a pd dataframe containing reads from all samples in one single chromosome 
    output: 
        coverage_result (dict): with key (sample) - value (average depth) pair
        empty_result (dict): with key(sample) - value(empty rate)
                """
    vcf_df = all_vcf_df.copy()
    coverage_result = {}
    empty_result = {}

    total_non_zero_count = 0
    zero_count = 0
    coverage_count = 0
    total_count = len(vcf_df.index) #total number of positions

    for sample in list(vcf_df.columns)[2:]:
        #check empty places
        vcf_df[sample] = vcf_df[sample].astype(str).str.replace('*','0', regex = True)
        zero_count = vcf_df[sample].astype(str).str.count('0').sum()
        empty_result[sample] = round(zero_count/total_count,3)

        total_non_zero_count = total_count - zero_count

        #check coverage = sum of length of reads / total_non_zero_count
        vcf_df[sample] = vcf_df[sample].astype(str).str.replace('0','', regex = True)
        coverage_count = vcf_df[sample].str.len().sum()
        coverage_result[sample] = round(coverage_count/total_non_zero_count,3)

    print('coverage and empty result done!')

    return coverage_result,empty_result

def check_nucleotide_diffs(all_vcf_df):
    """
    This function is used to check nucleotide identity difference between two mapping approaches
    Input:
        all_vcf_df: a pd dataframe containing reads from all samples in one single chromosome 
        
    """
    vcf_df = all_vcf_df.copy()
    print('copy done!')
    total_pos = len(vcf_df.index)
    mammoths = ['E467','L163','M6','oimyakon'] #specify sample (mammoth id
    compare_list = ['bwa-aln_vs_vg-all','bwa-aln_vs_vg-mammoth','vg-all_vs_vg-mammoth'] #specify comparision sequences
    diff_count= {}
    diff_result = {}

    for i in range(len(mammoths)): #which mammoth
        #compare each mappings result againt each other
        #convert reads to set then compare against eachother
        #generate a boolen list: equal - True, not equal - False

        #bwa-aln vs vg-all
        vcf_df[compare_list[0] + '_' + str(mammoths[i])] = np.where((vcf_df.iloc[:,i+2].map(set) == vcf_df.iloc[:,i+2+4].map(set)), True, False)
        diff_count[compare_list[0] + '_' + str(mammoths[i])] = (~vcf_df[compare_list[0] + '_' + str(mammoths[i])]).sum()

        #bwa-aln vs vg-mammoth
        vcf_df[compare_list[1] + '_' + str(mammoths[i])] = np.where((vcf_df.iloc[:,i+2].map(set) == vcf_df.iloc[:,i+2+4+4].map(set)), True, False)
        diff_count[compare_list[1] + '_' + str(mammoths[i])] = (~vcf_df[compare_list[1] + '_' + str(mammoths[i])]).sum()

        #vg-all vs vg-mammoth
        vcf_df[compare_list[2] + '_' + str(mammoths[i])] = np.where((vcf_df.iloc[:,i+2+4].map(set) == vcf_df.iloc[:,i+2+4+4].map(set)), True, False)
        diff_count[compare_list[2] + '_' + str(mammoths[i])] = (~vcf_df[compare_list[2] + '_' + str(mammoths[i])]).sum()

        print('this mammoth {mammoth} is done!'.format(mammoth = mammoths[i]))

    for comparision, count in diff_count.items():
        diff_rate = round(count/total_pos,3)
        diff_result[comparision] = diff_rate
    return diff_result

def check_empty(all_vcf_df):
    vcf_df = all_vcf_df.copy()
    empty_result = {}
    row_num = len(vcf_df.index)

    for sample in list(vcf_df.columns)[2:]:
        vcf_df[sample] = vcf_df[sample].astype(str).str.replace('*','0', regex = True)
        vcf_df[sample] = vcf_df[sample].astype(str).str.count('0')
        empty_result[sample] = round(vcf_df[sample].sum()/row_num,3)
    print('empty result done!')
    return empty_result

def report_result(o_file, result):
    with open(o_file, 'w') as f:
        for sample, coverage in result.items():
            f.write(str(sample) + '\t' + str(coverage) + '\n')
    f.close()



all_vcf_df = read_file_to_df(mpfile_file,index_file)

coverage_result,empty_result = check_coverage_empty_rate(all_vcf_df)
report_result('{chr}_coverage.txt'.format(chr = id),coverage_result)
report_result('{chr}_emptyrate.txt'.format(chr = id), empty_result)

diff_result = check_nucleotide_diffs(all_vcf_df)
report_result('{chr}_nuc_diff.txt'.format(chr = id), diff_result)

print('{chr} all check is done!'.format(chr = id))
