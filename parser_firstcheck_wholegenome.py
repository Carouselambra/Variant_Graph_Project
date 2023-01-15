from __future__ import division
import sys
import pandas as pd
import numpy as np
import gzip

"""This script is used to check coverage, nucleotide differences and empty mappings rates in the whole genome

Input:
    mpfile_file path
    index_file path

Output:
    three files:
        whole_genome_coverage
        whole_genome_emptyrate
        whole_genome_nuc_diff

Functions:
    read_file(mpfile, index file):
        open mpfile as a dataframe with columns name provided by index file
    report_result(output file name, result):
        write result (dict) to a txt file
    check_whole_coverage(index_file, chr_list):
        check the whole genome coverage of each sample (adding all chromosomes' coverages / non empty positions)
    check_whole_emptyrate(index_file, chr_list):
        check the whole genmoe empty rate of each sample (adding all chromosomes' empty positions /  total postions)
    check_whole_nuc_diff(index_file, chr_list):
        check the whole genome nucleotide differences between different mapping strategies (adding all chromosomes nucleotide differences between samples / total positions)

"""
#read index file
index_file = sys.argv[1]

#generate a list of chromosome 
chr_list = []
for i in range(1,28):
    chr_list.append('chr'+str(i))
chr_list.append('chrX')

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
    return all_vcf_df

def report_result(o_file, result):
    """this function is to convert a dictionary to a txt file
    input:
        o_file (str): output file path
        result (dict)

    output:
        a file with o_file path"""
    with open(o_file, 'w') as f:
        for sample, values in result.items():
            f.write(str(sample) + '\t' + str(values) + '\n')
    f.close()

def check_whole_depth(index_file, chr_list):
    """This function is used to check the average depth of all mapped samples
    input:
        index_file: filepath of the index file
        chr_list: a list of chromsome
    output:
        final_full_coverage_rate (dict) a dictionary with key (sample_name) - value (the corresponding average depth)
    """
    full_coverage = {} # [coverage reads count, non zero count] <= accumulate through each rounds
    final_full_coverage_rate = {} # report the final result

    #read index file
    with open(index_file,'r') as index:
        for line in index:
            line = line.strip()
            splited = line.split('\t')
            full_coverage[splited[1]+'_'+splited[2]] = []

    #read through each mpile file by chromosome
    for chr in chr_list:
        mpfile = "merged.chr{id}.mpile.gz".format(id = chr)
        this_chr_df = read_file(mpfile, index_file)
        total_count = len(this_chr_df.index)
        for sample in full_coverage:
            zero_count = 0
            non_zero_count = 0
            this_chr_df[sample] = this_chr_df[sample].astype(str).str.replace('*','0', regex = True)
            zero_count = this_chr_df[sample].astype(str).str.count('0').sum()
            non_zero_count = total_count - zero_count

            this_chr_df[sample] = this_chr_df[sample].astype(str).str.replace('0','', regex = True)
            coverage_count = this_chr_df[sample].str.len().sum()

            if full_coverage[sample] == []:
                full_coverage[sample].append(coverage_count)
                full_coverage[sample].append(non_zero_count)
            else:
                full_coverage[sample][0] += coverage_count
                full_coverage[sample][1] += non_zero_count
        print("this {sample} on {chromosome} coverage is added!".format(sample = sample, chromosome = chr))

    for sample, values in full_coverage.items():
        final_full_coverage_rate[sample] = round(values[0]/values[1],3)
    return final_full_coverage_rate

def check_whole_depth_distribution(index_file,chr_list):
    """this function is to show the overall depth distribution
    input: 
        index_file
            file path to the index file
        chr_list
            a list of chromosomes 
        bins 
            specify numbers of groups of 
    output:
        depth_distr (nested dict)
            a dictionary contains occurances of each depth of all samples
        std_dict 
            a dictionary contains standard deviations of all samples"""
    depth_distr = {}

     with open(index_file,'r') as id:
        sample_name = ""
        for line in id:
            line = line.strip()
            sample_name = line.split("\t")[1] + "_" + line.split("\t")[2]
            depth_distr[sample_name]

    

def check_whole_emptyrate(index_file, chr_list):
    full_emptyrcount = {}
    final_full_emptyrate = {}
    with open(index_file,'r') as index:
        for line in index:
            line = line.strip()
            splited = line.split('\t')
            full_emptyrcount[splited[1]+'_'+splited[2]] = []

    for chr in chr_list:
        mpfile = "merged.chr{id}.mpile.gz".format(id = chr)
        this_chr_df = read_file(mpfile, index_file)
        total_count = len(this_chr_df.index)

        for sample in full_emptyrcount:
            this_chr_df[sample] = this_chr_df[sample].astype(str).str.replace('*','0', regex = True)
            zero_count = this_chr_df[sample].astype(str).str.count('0').sum()

            if full_emptyrcount[sample] == []:
                full_emptyrcount[sample].append(zero_count)
                full_emptyrcount[sample].append(total_count)
            else:
                full_emptyrcount[sample][0] += zero_count
                full_emptyrcount[sample][1] += total_count
    print("this {sample} on {chromosome} empty count is added!".format(sample = sample, chromosome = chr))

    for sample, values in full_emptyrcount.items():
        final_full_emptyrate[sample] = round(values[0]/values[1],3)
    return final_full_emptyrate

def check_whole_nuc_diff(index_file, chr_list):
    full_nuc_diff_count = {}
    final_full_nuc_diff_rate = {}
    mammoths = ['E467','L163','M6','oimyakon'] #specify sample (mammoth id)
    compare_list = ['bwa-aln_vs_vg-all','bwa-aln_vs_vg-mammoth','vg-all_vs_vg-mammoth']
    for comparison in compare_list:
        for mammoth in mammoths:
            sample = comparison + '_' + mammoth
            full_nuc_diff_count[sample] = []

    for chr in chr_list:
        mpfile = "merged.chr{id}.mpile.gz".format(id = chr)
        this_chr_df = read_file(mpfile, index_file)
        total_pos = len(this_chr_df.index)
        diff_count = {} #temporary store for each chromosome

        for i in range(len(mammoths)):#which mammoth
            #compare each mappings result againt each other
            #convert reads to set then compare against eachother
            #generate a boolen list: equal - True, not equal - False

            this_chr_df[compare_list[0] + '_' + str(mammoths[i])] = np.where((this_chr_df.iloc[:,i+2].map(set) == this_chr_df.iloc[:,i+2+4].map(set)), True, False)
            diff_count[compare_list[0] + '_' + str(mammoths[i])] = (~this_chr_df[compare_list[0] + '_' + str(mammoths[i])]).sum()

            #bwa-aln vs vg-mammoth
            this_chr_df[compare_list[1] + '_' + str(mammoths[i])] = np.where((this_chr_df.iloc[:,i+2].map(set) == this_chr_df.iloc[:,i+2+4+4].map(set)), True, False)
            diff_count[compare_list[1] + '_' + str(mammoths[i])] = (~this_chr_df[compare_list[1] + '_' + str(mammoths[i])]).sum()

            #vg-all vs vg-mammoth
            this_chr_df[compare_list[2] + '_' + str(mammoths[i])] = np.where((this_chr_df.iloc[:,i+2+4].map(set) == this_chr_df.iloc[:,i+2+4+4].map(set)), True, False)
            diff_count[compare_list[2] + '_' + str(mammoths[i])] = (~this_chr_df[compare_list[2] + '_' + str(mammoths[i])]).sum()

        for combo, value in diff_count.items():
            if combo in full_nuc_diff_count:
                if full_nuc_diff_count[combo] == []:
                    full_nuc_diff_count[combo].append(value)
                    full_nuc_diff_count[combo].append(total_pos)
                else:
                    full_nuc_diff_count[combo][0] += value
                    full_nuc_diff_count[combo][1] += total_pos
        print("all mammoths has been compared on {chromosome}".format(chromosome = chr))

    for sample, values in full_nuc_diff_count.items():
        final_full_nuc_diff_rate[sample] = round(values[0]/values[1],3)

    return final_full_nuc_diff_rate



def check_wg_nuc_diff_exclude_empty(index_file,chr_list):

    """compare if nucleotide identity is same between two mapping approaches line by line
        note: both mapping approaches must have reads"""

    nuc_diff_dict = {} #key: mammmoths name, values: a dictionary with comparison # as key,  value is a list with non empty occurance and difference occurance
    mammoth_namelist = []
    #sample_names = []
    compare_mapping_index = ['bwa-aln_vs_vg-all','bwa-aln_vs_vg-mammoth','vg-all_vs_vg-mammoth']


    with open(index_file,'r') as id:
        mammoth_name = ""
        #sample_name = ""
        for line in id:
            line = line.strip()
            mammoth_name = line.split("\t")[2]
            if (mammoth_name in mammoth_namelist) == False:
                mammoth_namelist.append(mammoth_name)


    for i in mammoth_namelist:
        nuc_diff_dict[mammoth_namelist.index(i)] = {0:[0,0],1:[0,0],2:[0,0]}
    print('index copying done!')

    for chromosome in chr_list:
        mpfile = "merged.{chrname}.mpile.gz".format(chrname = chromosome)
        with gzip.open(mpfile, "r") as f1:
            for line in f1:
                line = line.decode()
                line = line.strip()
                splited = line.split("\t")

                for i in range(2,2+len(mammoth_namelist)): #iterate by mammoth

                    #bwa-aln vs vg_all
                    #ensure both mapping have reads
                    if (splited[i] != '*') and (splited[i+4] != '*'): 

                        #counting for shared read sites between bwa-aln and vg-all (comparison mapping index 0)
                        nuc_diff_dict[i-2][0][0] += 1

                        #check if both mapping have different nucleotides 
                        if set(splited[i]) != set(splited[i+4]):

                            #counting for read sites with different nucleotides
                            nuc_diff_dict[i-2][0][1] += 1 

                    #bwa-aln vs vg-mm
                    #ensure both mapping have reads
                    if (splited[i] != '*') and (splited[i+8] != '*'): 

                        #counting for shared read sites between bwa-aln and vg-mm (comparison mapping index 1)
                        nuc_diff_dict[i-2][1][0] += 1 

                        #check if both mapping have different nucleotides
                        if set(splited[i]) != set(splited[i+8]): 

                            #counting for read sites with different nucleotides
                            nuc_diff_dict[i-2][1][1] += 1 


                    #vg_all vs vg_mm
                    #ensure both mapping have reads
                    if (splited[i+4] != '*') and (splited[i+8] != '*'): 

                        #counting for shared read sites between vg_all and vg-mm (comparison mapping index 2)
                        nuc_diff_dict[i-2][2][0] += 1 

                        #check if both mapping have different nucleotides
                        if set(splited[i+4]) != set(splited[i+8]): 

                            #counting for read sites with different nucleotides
                            nuc_diff_dict[i-2][2][1] += 1 
        print(chromosome + 'is done!')
    #calculate the differences rates and store them into a new result dictionary
    entry_name = ""
    diff_rate = 0
    nuc_diff_rate_dict = {}
    for i in range(4):
        for j in range(3):
            entry_name = ""
            diff_rate = 0
            entry_name = compare_mapping_index[j] + '_' + mammoth_namelist[i]
            diff_rate = nuc_diff_dict[i][j][1] / nuc_diff_dict[i][j][0]
            nuc_diff_rate_dict[entry_name] = diff_rate

    return nuc_diff_rate_dict


    
#wg_coverage_result = check_whole_depth(index_file, chr_list)
#report_result("whole_genome_coverage", wg_coverage_result)
#print("whole_genome_coverage done!")

#wg_emptyrate_result = check_whole_emptyrate(index_file,chr_list)
#report_result("whole_genome_emptyrate", wg_emptyrate_result)
#print("whole_genome_emptyrate done!")

#wg_nuc_diff_rate_result = check_whole_nuc_diff(index_file, chr_list)
#report_result("whole_genome_nuc_diff",wg_nuc_diff_rate_result)
#print("whole_genome_nuc_diff done!")

wg_nuc_diff_rate_exclude_0 = check_wg_nuc_diff_exclude_empty(index_file,chr_list)
report_result("whole_genome_nuc_diff_exclude_empty", wg_nuc_diff_rate_exclude_0)
print('whole_genome_nuc_diff_exclude_empty done!')
