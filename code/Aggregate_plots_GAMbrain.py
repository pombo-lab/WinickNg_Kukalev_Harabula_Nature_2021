get_ipython().run_line_magic('matplotlib', 'inline')
import pandas as pd
import random
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import scipy as spagrg
from scipy import stats
from scipy import ndimage
from matplotlib import transforms
import os, sys
import glob
import seaborn as sb

"""
Aggregate_plots_GAMbrain.py
Script written by: Warren Winick-Ng

This script is part of the study:
Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded by specific chromatin topologies", Nature 2021 

Date: August 16, 2021


This script takes an input set of contacts, and will generate 9x9 bin aggregate plots of the contact and surrounding region.

The main function is 'get_aggregate_plot_forfile':
'path_for_input': Each .bed input file in folder loc should be in long format: 'chr1', 'start1', 'stop1', 'chr2', 'start2', 'stop2', 'score' *note: remove 'score' if using randomized contacts*; 
'NPMI_folder': specify the folder for the NPMI matrices for the dataset of interest;
'outpath': specify the folder to save the aggregate maps'''

This code was generated in and works best if run in the Jupyter environment

Usage example:
>>> get_aggregate_plot_forfile('path_for_input_dataframe', 'path_for_NPMI_values_for_each_chromosome', 'path_to_output')

Please refer to the Read.me document for details on usage and distribution. 

"""

def get_resolution (matrix):
    index_bin = matrix.columns[0]
    coord = index_bin.split(':')[1]
    start, stop = coord.split ('-')
    resolution = int (stop) - int (start)
    return resolution

def get_chromosome_matrix (matrices_folder, chromosome):
    matrices_list = glob.glob(matrices_folder + '*.txt.gz')
    for matrix_file in matrices_list:
        file_name = matrix_file.split('/')[-1]
        file_name_split = file_name.split('.')
        if chromosome in file_name_split:
            print ('Working on chromosome: ' + chromosome + ',', end=" ")
            matrix = pd.read_csv (matrix_file, sep='\t', index_col=[0])
        if not matrix_file:
            sys.exit('Can not find specified chromosome in the folder')
    return matrix

def zscore_matrix (NPMI_matrix):
    """
    Calculate Z-score for each diagonal of an input matrix
    :returns: normalized z-score matrix with index and column names from orginal file
    """
    chrom_matrix = np.ma.masked_invalid(NPMI_matrix)
    z_matrix = np.zeros_like(chrom_matrix)
    for i in range(1, len(chrom_matrix)-1):
        diag = np.diagonal(chrom_matrix, offset=i) 
        zdiag = stats.zscore(diag) 
        np.fill_diagonal(z_matrix[:,i:], zdiag)
    Z_matrix = z_matrix + z_matrix.T
    return pd.DataFrame(Z_matrix, index=NPMI_matrix.columns, columns=NPMI_matrix.columns)

def aggregated_maps (df_longbed, matrices_folder, delta_bins=4):
    sumup_matrix_forfinal = np.zeros ([2*delta_bins + 1, 2*delta_bins + 1]).astype(int) # Empty matrix of expected shape
    chromosomes = list (set (df_longbed ['chrom1'].tolist()))
    chromosomes.sort()
    for chromosome in chromosomes:
        sumup_matrix = np.zeros ([2*delta_bins + 1, 2*delta_bins + 1]).astype(int)
        subset_contacts_df = df_longbed.query ('chrom1 == @chromosome')
        entire_chromosome_matrix = get_chromosome_matrix(matrices_folder, chromosome)
        entire_chromosome_Zmatrix = zscore_matrix(entire_chromosome_matrix)
        resolution = get_resolution (entire_chromosome_matrix)
        chromosome_size = len (entire_chromosome_matrix)
        for index, row in subset_contacts_df.iterrows():
            chrom, start1, stop1, start2, stop2  = row[0], row[1], row[2], row[4], row[5]
            start_bin1 = int (start1 / resolution) - delta_bins
            stop_bin1 = int (stop1 / resolution) + delta_bins
            start_bin2 = int (start2 / resolution) - delta_bins
            stop_bin2 = int (stop2 / resolution) + delta_bins
            if stop_bin2 < chromosome_size and start_bin1 > 0:# Check that the region is not at the beginning or end of the chromosome
                if (start_bin2 - start_bin1) > delta_bins: # The interacting windows should be at-least 3 bins away 
                    matrix_subset = entire_chromosome_Zmatrix.iloc [start_bin1:stop_bin1, start_bin2:stop_bin2]
                    #boolean_matrix = matrix_subset.notna() #use this if you want a count matrix instead of values
                    values_matrix = matrix_subset.values
                    values_matrix[np.isnan(values_matrix)] = 0
                    sumup_matrix = np.add(sumup_matrix, values_matrix) 
        sumup_matrix_mean_chrom = sumup_matrix / (int(len(subset_contacts_df)))
        sumup_matrix_mean_chrom[np.isnan(sumup_matrix_mean_chrom)] = 0
        sumup_matrix_forfinal = np.add(sumup_matrix_forfinal, sumup_matrix_mean_chrom)
    sumup_matrix_mean = sumup_matrix_forfinal / (int(len(chromosomes)))
    return sumup_matrix_mean



def get_aggregate_plot_forfile(path_for_input, NPMI_folder, outpath):
    '''Aggregate plotting script for NPMI input data from a folder; \\
    'path_for_input': Each .bed input file in folder loc should be in long format: \\
    'chr1', 'start1', 'stop1', 'chr2', 'start2', 'stop2', 'score' *NOTE: remove 'score' if using randomized contacts*; \\
    'NPMI_folder': specify the folder for the NPMI matrices for the dataset of interest; \\
    'outpath': specify the folder to save the aggregate maps'''
    
    dirs = os.listdir(path_for_input)
    file_list = []

    for file in dirs:
        file_list.append(file)
    
    for file in file_list:
        file_df = pd.read_csv(path_for_input + file, sep='\t')
        file_name = str(file)
        file_name_short = file_name[:-4] #format for .bed files, modify if different file type
        file_df.columns = ['chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2', 'score'] #remove 'score' form this line if not in file
        file_df_input = file_df[['chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2']]
        
        print('Starting to work on aggregate map of: ' + file_name)        
        aggre_matrix = aggregated_maps (file_df_input, NPMI_folder)
        print('finished aggregate map of: ' + file_name)
        
        #plot figure
        plt.figure (figsize=(12,12))
        matrix1_plt = plt.imshow (aggre_matrix, interpolation=None, cmap='RdYlBu_r', vmin=0, vmax=1)
        plt.colorbar(matrix1_plt, orientation='horizontal', shrink = 0.5, aspect=15)
        plt.tight_layout()
        plt.savefig(outpath + file_name_short + '_aggregatemap.png', dpi=300)
        plt.close()
        print('Saved aggregate map of: ' + file_name + '\n')
    print('All files saved!')

