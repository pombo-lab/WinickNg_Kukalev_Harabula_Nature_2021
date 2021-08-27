#!~/.guix-profile/bin/python3

"""
====================================================================================================================
Python script name GAM.NPMI.coseg.matrices.py
Script written by Alexander Kukalev

This script is part of the study: Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded 
by specific chromatin topologies", Nature 2021

Date: August 16, 2021

This script produce raw co-segregation or NPMI normalized pair-wise contact matrices for all segregation table in the 
specified input folder. NPMI matrices can also be Z-score normalized. The matrices for each chromosome are saved as separate files

Usage example:
>>> python3 GAM.NPMI.coseg.matrices.py $input_folder $output_folder $matrixtype 

$matrixtype can take three values 'NPMI', 'coseg' or 'Z' for normalized NPMI, 
raw co-segregation or Z-score normalized NPMI matrices respectively

=====================================================================================================================
"""
#!~/.guix-profile/bin/python3

import sys
import pandas as pd
import numpy as np
from scipy import stats

def get_chromosomes_list (segregation_table):
    # This function extract the list of chromosomes from GAM segregation table
    # Random chromosomes and mitotic chromosome are excluded
    all_chromosomes = list (set (segregation_table.index.get_level_values('chrom')))
    all_chromosomes.sort()
    karyotype_set = []
    for chromosome in all_chromosomes:
        if (chromosome.find ('_') == -1) and (chromosome != "chrM"):
            karyotype_set.append (chromosome)
    return karyotype_set

def calculate_NPMI(seg_matrix):
    # This function calculate NPMI normalized pair-wise chromatin contact matrix from GAM segregation table
    np.warnings.filterwarnings('ignore')
    # calculate M the number of samples
    M = len(seg_matrix[0])
    #calculate p(x,y) matrix i.e coseg matrix divided my M
    pxy = (seg_matrix.dot(seg_matrix.T))/M
    #p(x)p(y) matrix  - i.e detection frequency with the dot product of it's transposed self equals an N*N matrix
    pxpy = seg_matrix.sum(1).reshape(-1,1)/M * seg_matrix.sum(1)/M
    #define PMI matrix
    PMI = np.log2(pxy/pxpy)
    #bound values between -1 and 1
    NPMI = PMI/-np.log2(pxy)
    return NPMI

def zscore_matrix (NPMI_matrix):
    # This function calculate Z-score normalized matrix from NPMI normalized pair-wise matrix
    chrom_matrix = NPMI_matrix.copy()
    z_matrix = np.zeros_like(chrom_matrix)
    for i in range(1, len(chrom_matrix)-1):
        diag = np.diagonal(chrom_matrix, offset=i) 
        zdiag = stats.zscore(diag, nan_policy='omit') 
        np.fill_diagonal(z_matrix[:,i:], zdiag)
    Z_matrix = z_matrix + z_matrix.T
    return Z_matrix

def calculate_cosegregation(seg_matrix):
    # This function calculate raw co-segregation matrices from GAM segregation table
    M = len(seg_matrix[0])
    pxy = (seg_matrix.dot(seg_matrix.T))/M
    return pxy

def generate_matrices_for_each_segregation_table_in_a_folder (input_folder, output_folder, matrix_type):
    # This is the main function 
    import glob
    seg_tables_list = glob.glob(input_folder + '*.table')
    for seg_table in seg_tables_list:
        print ('Reading segregation table ' + seg_table)
        seg_table_df = pd.read_csv (seg_table, sep='\t', index_col=[0,1,2])
        seg_table_name = seg_table.split('/')[-1]
        seg_table_name = seg_table_name.replace ('.table', '')
        list_of_chromosomes = get_chromosomes_list (seg_table_df) 
        for chrom in list_of_chromosomes:
            print ('Generating matrices for chromosome ' + chrom)
            subset_segregation = seg_table_df.loc[chrom]
            seg_matrix = subset_segregation.values
            if matrix_type=='NPMI':
                matrix = calculate_NPMI (seg_matrix)
                matrix_file_name = output_folder + str(seg_table_name) + '.NPMI.' + chrom + '.txt.gz'
            if matrix_type=='coseg':
                matrix = calculate_cosegregation (seg_matrix)
                matrix_file_name = output_folder + str(seg_table_name) + '.coseg.' + chrom + '.txt.gz'
            if matrix_type=='Z':
                matrix_file_name = output_folder + str(seg_table_name) + '.NPMI.' + chrom + '.txt.gz'
                matrix = calculate_NPMI (seg_matrix)
                Z_matrix = zscore_matrix (matrix)
            matrix_df = pd.DataFrame(matrix,index=subset_segregation.index, columns=subset_segregation.index) 
            # Get coordinates for index and header
            subset_segregation.reset_index(inplace = True)
            coord = chrom + ':' + subset_segregation['start'].astype(str) + '-' + subset_segregation['stop'].astype(str)
            matrix_df = pd.DataFrame(matrix,index=coord, columns=coord) 
            print ('Saving ' + matrix_type + ' matrix as ' + matrix_file_name)
            matrix_df.to_csv (matrix_file_name, sep='\t', compression='gzip', na_rep='NA')
            if matrix_type=='Z':
                Z_score_df = pd.DataFrame(Z_matrix,index=coord, columns=coord) 
                Z_matrix_file_name = output_folder + str(seg_table_name) + '.Z_score.' + chrom + '.txt.gz'
                print ('Saving Z-score matrix as ' + Z_matrix_file_name)
                Z_score_df.to_csv (Z_matrix_file_name, sep='\t', compression='gzip', na_rep='NA')

#========================================================================================

if len(sys.argv)==4:
    (input_folder, output_folder, matrix_type) = sys.argv[1:]
else:
    print ('You need to provide the path to the input folder, path to the output_folder and type of the matrices (npmi or coseg)')
    sys.exit()

generate_matrices_for_each_segregation_table_in_a_folder (input_folder, output_folder, matrix_type)


