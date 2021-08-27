#!~/.guix-profile/bin/python3

"""
====================================================================================================================
Python script name GAM.segregation.tables.and.stats.py
Script written by Alexander Kukalev

This script is part of the study: Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded 
by specific chromatin topologies", Nature 2021

Date: August 16, 2021

This script generate GAM segregation tables from the coverage tables and produce the GAM QC tables with most important 
metrics for each sample. Script GAM.segregation.tables.and.stats.sh will generate GAM reads coverage tables for multiple
resolutions and automatically execute this script

Usage example: 
see script GAM.segregation.tables.and.stats.sh 
>>> python3 GAM.segregation.tables.and.stats.py $FOLDER $DATASET_NAME $COVERAGE_TABLE $RESOLUTION   

By default, the script will  produce QC stats for every sample at 50Kb resolution
=====================================================================================================================
"""

import sys
import numpy as np
import pandas as pd
from datetime import date

""" All functions """ 

""" Get mapping stats from fastq and bam files """ 
def get_mapping_stats (import_folder):
    # This function will calculate reads dependent QC metrics for each sample
    # such as total number of reads, number of mapped reads, number of uniquely mapped reads,
    # and percent of mapped reads to the target genome
    print ('Starting to calculate mapping stats...')
    import glob
    import gzip
    import pysam
    fastqlist = glob.glob(import_folder + '*.fastq.gz')
    bamlist = glob.glob(import_folder + '*.sorted.bam')
    rmduplist  = glob.glob(import_folder + '*.rmdup.bam')
    fastqlist.sort()
    bamlist.sort()
    rmduplist.sort()
    samples_list, total_reads, mapped_reads, unique_reads, percent_list = [], [], [], [], []
    for fastq_file, bam_file, rmdup_file in zip (fastqlist, bamlist, rmduplist):
        samples_list.append (fastq_file)
        bamf = pysam.Samfile(bam_file, 'rb')
        rmdupf = pysam.Samfile(rmdup_file, 'rb')
        bam_reads = 0
        number_of_reads = 0
        rmdup_reads = 0
        with gzip.open(fastq_file, 'rb') as fastq:
            for id in fastq:
                seq = next(fastq)
                number_of_reads += 1
                next(fastq)
                next(fastq)
        for read in bamf.fetch():
            bam_reads += 1
        for read in rmdupf.fetch():
            rmdup_reads += 1
        percent = (float(bam_reads/number_of_reads))*100
        total_reads.append (number_of_reads)
        mapped_reads.append (bam_reads)
        unique_reads.append (rmdup_reads)
        percent_list.append (percent)
    stats_df = pd.DataFrame(list(zip(fastqlist, rmduplist, total_reads, mapped_reads, unique_reads, percent_list)), \
                            columns=['Fastq_file',  'Rmdup_path', 'Total reads', 'Mapped_reads', 'Uniquely_mapped_reads', 'Percent_of_mapped_reads'])
    print ('Mapping stats are done')
    return stats_df

""" Generate segregation table from coverage table"""
def remove_path_from_columns_names (df):
    # This function simplify column names in coverage tables by removing full path
    list_of_columns = df.columns.tolist()
    simplified_columns = []
    for file_name in list_of_columns:
        if file_name.find('/') != -1:  # Remove path from the file name is present
            file_name = file_name.split('/')[-1]
            file_name = file_name.split ('.')[0]
        else:
            file_name = file_name.split ('.')[0]
        if file_name[0].isdigit(): # Add letter to the file name is starts from digit
            file_name = 'f_' + file_name
        if file_name.find('-') != -1: # Replace '-' to '_' if present
            file_name = file_name.replace ('-', '_')
        simplified_columns.append (file_name)
    return_df = df.copy()
    return_df.columns = simplified_columns
    return return_df

def get_resolution (segregation_table):
    # Identify resolution of coverage or segregation table
    start_values = np.array (segregation_table.index.get_level_values('start'))
    resolution = start_values [1] - start_values [0] # Find the resolution from table 
    return resolution 

def get_lowest_percentile (coverage_df):
    # Defined based on the resolution
    # For resolutions below 250Kb the algorithm scan all distribution from 0 to 98 percentile
    # For resolutions above 250Kb  the algorithm scan distribution from 40 to 98 percentile
    resolution = get_resolution (coverage_df)
    if resolution >= 250000: # Set percentile boundaries depending on the resolution
        lowest_percentile = 40 # For resolutions above 250Kb go from 40 to 98 percentile
    else:
        lowest_percentile = 0 # Fr resolutions below 250Kb go from 0 to 98 percentile
    return lowest_percentile

def get_chromosomes_list (coverage_df):
    chromosomes_list = []
    chromosome_values = np.array (coverage_df.index.get_level_values('chrom'))
    chromosomes = set (chromosome_values)
    for chromosome in chromosomes:
        if (chromosome.find ('_') == -1) and (chromosome != 'chrM'):
            chromosomes_list.append (chromosome)
    chromosomes_list.sort()
    return chromosomes_list
     
def get_segregation (coverage_per_bin, threshold):
    segregation_per_bin = coverage_per_bin > threshold
    segregation_per_bin = segregation_per_bin.astype(int)
    return segregation_per_bin
             
def get_orphan_windows (segregation_per_bin, chromosomes):
    orphan_windows = []
    aux, counter = 0, 0
    # count orphan windows
    for chromosome in chromosomes:
        value_list = segregation_per_bin.loc[chromosome].values
        for i,window in enumerate(value_list):
            if window == 1:
                counter += 1
                if i == 0:
                    if sum(value_list[i:i+2]) == 1:
                        aux += 1
                elif i == len(value_list):
                    if sum(value_list[i-1:i+1]) == 1:
                        aux += 1
                else:
                    if sum(value_list[i-1:i+2]) == 1:
                        aux += 1     
    try:
        orphan_percentage = (aux/float(counter))*100
    except:
        orphan_percentage = 0
    return orphan_percentage

def read_coverage_table (coverage_table_path):
    print ('Reading coverage table...')
    coverage_df = pd.read_csv (coverage_table_path, sep='\t', index_col=[0,1,2])
    coverage_df_simple_names = remove_path_from_columns_names (coverage_df)
    return coverage_df_simple_names

def get_segregation_table (coverage_table_path, dataset_name):
    coverage_df = read_coverage_table (coverage_table_path)
    chromosomes = get_chromosomes_list (coverage_df) # Get the list of chromosomes
    lowest_percentile = get_lowest_percentile (coverage_df) # Get the lowest percentile depending on the resolution
    list_of_samples = coverage_df.columns.tolist() 
    segregation_df = pd.DataFrame(index=coverage_df.index)
    print ('Starting to generate segregation table...\n')
    for sample in list_of_samples:
        #sys.stdout.write(".")
        #sys.stdout.flush()
        coverage_per_bin = coverage_df [sample]
        coverage_above0 = coverage_per_bin[coverage_per_bin > 0].tolist()
        coverage_per_bin_log10 = np.log10(coverage_above0)
        list_of_orphans, list_of_percentiles = [], []
        if coverage_per_bin_log10.size != 0: # Calculate optimal threshoold for samples that have reads and coverage 
            biggest_percentile = 98
            for percentile in reversed(range(lowest_percentile, biggest_percentile + 1)):
                threshold = np.percentile (coverage_per_bin_log10, percentile)
                threshold_in_nucleotides = pow (10, threshold)
                sample_segregation = get_segregation (coverage_per_bin, threshold_in_nucleotides)
                percent_of_orphan_windows = get_orphan_windows (sample_segregation, chromosomes)
                list_of_orphans.append (percent_of_orphan_windows)
            for pos in range (len(list_of_orphans)):
                if list_of_orphans[pos]==0:
                    list_of_orphans[pos]=100
            percentile_with_lowest_orphans = np.argmin(list_of_orphans)
            optimal_threshold = np.percentile (coverage_per_bin_log10, (biggest_percentile - percentile_with_lowest_orphans))
            optimal_threshold_estimate = pow (10, optimal_threshold)
            if optimal_threshold_estimate > 78:
                optimal_threshold_in_nucleotides = optimal_threshold_estimate
            else:
                optimal_threshold_in_nucleotides = 78
            # Calculate segregation per sample
            segregation_per_sample = get_segregation (coverage_per_bin, optimal_threshold_in_nucleotides)
            print ('Threshold for sample ' + sample + ' ' + str(int(optimal_threshold_in_nucleotides)) + ' bp')
        else: # For samples that do not have reads and coverage generate an empty column with zeros for all windows
            segregation_per_sample = pd.Series([0] * len(coverage_df.index))
        segregation_df [sample] = pd.DataFrame (segregation_per_sample, index=segregation_df.index)
    resolution = get_resolution (segregation_df)
    name_of_segregation_table = dataset_name + '.segregation_at' + str(resolution) + '.' + str(date.today()) + '.table'
    segregation_df.to_csv (name_of_segregation_table, sep='\t', index=True, header=True)
    print ('Segregation table was saved as ' + name_of_segregation_table)
    return segregation_df

"""Calculate total genome coverage per sample"""
def get_genome_coverage_per_sample (segregation_table): 
    # Calculate total genome coverage for all samples in the segregation table
    genome_coverage = []
    list_of_samples = segregation_table.columns.tolist()
    for sample in list_of_samples:
        number_of_positive_windows = len(segregation_table.query (f'{sample}==1'))
        total_windows = len (segregation_table)
        genome_coverage_per_sample = (float(number_of_positive_windows / total_windows) * 100)
        genome_coverage.append (genome_coverage_per_sample)
    print ('Genome coverage for each sample was calculated')
    return list_of_samples, genome_coverage

"""Calculate orphan windows per sample"""
def get_orphan_windows_per_sample (segregation_table):
    # Calculate percent of orphan windows (positive windows without neightboors) 
    # for all samples in the segregation table
    orphan_windows = []
    chromosomes = get_chromosomes_list (segregation_table) # Get the list of chromosomes
    list_of_samples = segregation_table.columns.tolist() # Get the list of all columns
    for sample in list_of_samples: # Iterate over samples
        windows_per_sample = segregation_table [sample]
        aux, counter = 0, 0
        # count orphan windows
        for chromosome in chromosomes:
            value_list = windows_per_sample.loc[chromosome].values
            for i,window in enumerate(value_list):
                if window == 1:
                    counter += 1
                    if i == 0:
                        if sum(value_list[i:i+2]) == 1:
                            aux += 1
                    elif i == len(value_list):
                        if sum(value_list[i-1:i+1]) == 1:
                            aux += 1
                    else:
                        if sum(value_list[i-1:i+2]) == 1:
                            aux += 1     
        try:
            orphan_percentage = (aux/float(counter))*100
        except:
            orphan_percentage = 0
        orphan_windows.append (orphan_percentage)
    print ('Percent of orphan windows for each sample was calculated')
    return orphan_windows

"""Calculate average number of windows in 20xresolution region per sample"""
def get_average_windows_in_region (segregation_table, windows_in_a_region=20):
    chromosomes = get_chromosomes_list (segregation_table) # Get the list of chromosomes    
    number_of_chromosomes = len (chromosomes)
    list_of_samples = segregation_table.columns.tolist() 
    average_window_in_region_per_sample, average_mean_windows_per_sample = [], []
    for sample in list_of_samples:
        subset = segregation_table[sample]
        average_mean_windows_per_chromosome, number_of_windows_in_region = [], []
        for chromosome in chromosomes:
            i = 0
            mean_of_windows = []
            while (i < len(subset.loc[chromosome].values)):
                values = subset.loc[chromosome].values[i:i+windows_in_a_region]
                if sum(values) != 0:
                    number_of_windows_in_region.append(sum(values))
                    mean_of_windows.append(np.mean(values))
                i += windows_in_a_region
            average_mean_windows_per_chromosome.append (mean_of_windows)
        if number_of_windows_in_region != []:
            average_window_in_region_per_sample.append (np.mean(number_of_windows_in_region))
        else:
            average_window_in_region_per_sample.append (0)
        average_mean_windows_per_sample.append (average_mean_windows_per_chromosome)
    resolution = get_resolution (segregation_table)
    print ('Average number of windows in ' + str(20*resolution) + ' bp was calculated')
    return average_window_in_region_per_sample, average_mean_windows_per_sample

"""Calculate number of chromosomes per sample"""
def get_number_of_chromosomes_per_sample (mean_windows_in_region):
    number_of_chromosome_per_sample = []
    for sample in mean_windows_in_region:
        flat_list = [item for sublist in sample for item in sublist]
        if flat_list != []:
            cutoff = np.percentile (flat_list, 75)
            number_of_positive_chromosomes = 0
            for chromosome in sample:
                 for region in chromosome:
                    if region  > cutoff:
                        number_of_positive_chromosomes+=1
                        break
            number_of_chromosome_per_sample.append (number_of_positive_chromosomes)
        else:
            number_of_chromosome_per_sample.append (0)
    print ('Average number of chromosomes in each sample was calculated')
    return number_of_chromosome_per_sample

###############################################################################

import_folder = sys.argv[1]
dataset_name = sys.argv[2]
coverage_table_path = sys.argv[3]
stats = sys.argv[4]

print ('Processing GAM samples located in the folder ' + import_folder)
print ('Name of the dataset ' + dataset_name)
print ('Segregation table will be generated from ' + coverage_table_path)

segregation_df = get_segregation_table (coverage_table_path, dataset_name)

if stats=='50000':
    stats_df = get_mapping_stats (import_folder)
    samples_list, coverage_list = get_genome_coverage_per_sample (segregation_df)
    orphans_list = get_orphan_windows_per_sample (segregation_df)
    number_of_windows_in_region, average_windows_mean = get_average_windows_in_region (segregation_df)
    number_of_chromosomes = get_number_of_chromosomes_per_sample (average_windows_mean)
    resolution = get_resolution (segregation_df)
    stats_df ['Simple_name'] = pd.DataFrame (samples_list, index=stats_df.index)
    stats_df ['Genome_coverage'] = pd.DataFrame (coverage_list, index=stats_df.index)
    stats_df ['Percent_of_orphans'] = pd.DataFrame (orphans_list, index=stats_df.index)
    stats_df ['Average number of windows in ' + str(20*resolution)] = pd.DataFrame (number_of_windows_in_region, index=stats_df.index)
    stats_df ['Number_of_chromosomes'] = pd.DataFrame (number_of_chromosomes, index=stats_df.index)
    stats_df.to_csv (dataset_name + '.stats.at_' + str(resolution) + '.' + str(date.today()) + '.table', index=False, header=True, sep='\t')
    print ('Stats table was saved as ' + dataset_name + '.stats.at_' + str(resolution) + '.' + str(date.today()) + '.table')