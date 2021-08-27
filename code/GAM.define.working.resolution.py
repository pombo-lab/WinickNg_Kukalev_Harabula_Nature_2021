#!~/.guix-profile/bin/python3

"""
====================================================================================================================
Python script name GAM.define.working.resolution.py
Script written by Alexander Kukalev

This script is part of the study: Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded 
by specific chromatin topologies", Nature 2021

Date: August 16, 2021

This script tests the quality of genome sampling for every GAM segregation table in a folder
and inform the user if there are enough samples to work at specified resolution

Usage example:
>>> python3 GAM.define.working.resolution.py $input_folder $genomic_distance 

To test the quality of genome sampling for all distances set genomic_distance to zero.
=====================================================================================================================
"""
def calculate_cosegregation (segregation_table):
    # The function calculate cosegregation matrix
    pxy = (segregation_table.dot(segregation_table.T))
    return pxy

def calculate_kendalltau (y_normal, y_data):
    # The function calculate kendal rank tau coefficient
    import scipy.stats
    tau = scipy.stats.kendalltau (y_normal, y_data)[0]
    return tau

def get_resolution (segregation_table):
    # The function get the resolution of provided segtable
    import pandas as pd
    import numpy as np
    if len (segregation_table.index.names) == 1: # Check if supplied segregation table has multiindex
        segregation_table.set_index(['chrom','start','stop'], inplace = True) 
    start_values = np.array (segregation_table.index.get_level_values('start'))
    resolution = start_values [1] - start_values [0] # Find the resolution from table 
    return resolution  

def remove_never_detected_windows (seg_table):
    # This function removes never detected unmappable genomic windows
    import numpy as np
    import pandas as pd
    windows_detection_across_dataset = seg_table.sum(axis=1).values
    seg_table ['WinSum'] = pd.DataFrame (windows_detection_across_dataset, index=seg_table.index)
    segtable_without_not_detected_windows = seg_table.query ('WinSum != 0 ')
    segtable_without_not_detected_windows.drop (columns=['WinSum'], inplace=True)
    return segtable_without_not_detected_windows

def get_segvalues_for_certain_genomic_distance (segregation_table, distance_from_diagonal, region):
    # The function subset segtable for single chromosome (or any other region), remove never detected windows, 
    # calculate co-segregation matrix and return an array with co-segregation values at particular genomic distance
    import pandas as pd
    import numpy as np
    matrix_values = []
    resolution = get_resolution (segregation_table)
    subset_table = segregation_table.loc [region]
    subset_table_remove_non_detected_windows = remove_never_detected_windows (subset_table)
    segmatrix = calculate_cosegregation (subset_table_remove_non_detected_windows)
    if distance_from_diagonal != 0:
        bins_from_diagonal = int (distance_from_diagonal/resolution)
    else: # do for entire matrix 
        bins_from_diagonal = len (segmatrix)
    for diagonal in range (1, bins_from_diagonal):
        diagonal_values = (np.diag (segmatrix, diagonal))
        matrix_values = np.append (matrix_values, diagonal_values)
    return matrix_values  
    
def yeojohnson_transform (input_aray):
    from scipy.stats import yeojohnson
    transformed, lmbda = yeojohnson(input_aray)
    return transformed

def qq_plot(data):
    import numpy as np
    sample_size = len (data)
    mean = np.mean (data)
    std = np.std (data)
    normal = np.random.normal(loc = mean, scale = std, size=sample_size)
    np.random.shuffle(data)
    np.random.shuffle(normal)
    x = np.sort(normal[0:sample_size])
    y_normal = np.sort(normal[0:sample_size])
    y_data = np.sort(data[0:sample_size])
    return x, y_normal, y_data

def test_normality_of_coseg_values_distribution (seg_table, distance, region = 'chr10'):
    # The main function, by default the test is done for chr10.
    import scipy.stats
    import numpy as np
    from colorama import Fore, Style
    if len (seg_table.index.names) == 1: # Check if supplied segregation table has multiindex
            seg_table.set_index(['chrom','start','stop'], inplace = True)
    # Since never detected windows will never co-segregate with any other genomic window and
    # are removed from the segtable, it's not advisable to perform this test with less than 100 NPs
    if len (seg_table.columns) < 100: 
        print ('Warning! The size of your GAM dataset seems to be too small for this test')
    segvalues = get_segvalues_for_certain_genomic_distance (seg_table, distance, region)
    segvalues_transformed = yeojohnson_transform (segvalues)
    x, y_normal, y_data = qq_plot (segvalues_transformed)
    tau = calculate_kendalltau (y_normal, y_data)
    lenght = len (segvalues)
    values_above_1 = np.count_nonzero (segvalues > 1)
    percent = (values_above_1 / lenght) * 100
    if tau > 0.95:
        print (f'{Fore.GREEN}Quality of genome sampling is good for specified distance and resolution{Style.RESET_ALL}')
        print (f'{Fore.GREEN}Percent of never detected or detected only ones window pairs --> {Style.RESET_ALL}' + str (100 - percent))
        print (f'{Fore.GREEN}Tau --> {Style.RESET_ALL}' + str (tau))
    else:
        print (f'{Fore.RED}You need more samples to be able to work at this genomic distance and resolution{Style.RESET_ALL}')
        print (f'{Fore.RED}Percent of never detected or detected only ones window pairs --> {Style.RESET_ALL}' + str (100 - percent))
        print (f'{Fore.RED}Tau --> {Style.RESET_ALL}' + str (tau))

#==========================================

import_folder = sys.argv[1]
distance = sys.argv[2]

import glob
seg_tables_list = glob.glob(input_folder + '*.table')
seg_tables_list.sort()
for seg_table in seg_tables_list:
    seg_table_df = pd.read_csv (seg_table, sep='\t', index_col=[0,1,2])
    seg_table_name = seg_table.split('/')[-1]
    print (seg_table_name)
    test_normality_of_coseg_values_distribution (seg_table_df, distance)
    print ()