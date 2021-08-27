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
Randomcontacts_fromlist_GAMbrain.py
Script written by: Warren Winick-Ng

This script is part of the study:
Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded by specific chromatin topologies", Nature 2021 

Date: August 16, 2021


This script takes an input set of contacts, and will generate:
-3 sets of randomized contacts matched for chromosome and distance

This code was generated in and works best if run in the Jupyter environment.
Output from this code was intended for use with: 'Aggregate_plots_GAMbrain.py'

Please refer to the Read.me document for details on usage and distribution. 

"""



def get_random_contacts(longdf):
    chrome_sizes = {"chr1" : 195471971,"chr2" : 182113224,"chr3" : 160039680,"chr4" : 156508116, "chr5" : 151834684,"chr6" : 149736546,"chr7" : 145441459,"chr8" : 129401213, "chr9" : 124595110,"chr10" : 130694993, "chr11" : 122082543,"chr12" : 120129022, "chr13" : 120421639,"chr14" : 124902244,"chr15" : 104043685,"chr16" : 98207768, "chr17" : 94987271,"chr18" : 90702639,"chr19" : 61431566,"chrX": 171031299,"chrY" : 91744698} 
    rand_stop = []

    for row in longdf.itertuples(index=True, name='Pandas'):
        chrom  = row.chrom1
        end_chromosome = int(chrome_sizes[str(chrom)]) - int(longdf['distance'].max())
        rand_stop_number = random.randint(5000000,end_chromosome)
        rand_stop.append(rand_stop_number)
    
    return rand_stop

#get 3 sets of randomized contacts as single coordinates from each input file in a folder
def get_random_coords_from_file(path_for_input, outpath):
    dirs = os.listdir(path_for_input)
    file_list = []

    for file in dirs:
        file_list.append(file)
 
    for file in file_list:
        file_df = pd.read_csv(path_for_input + file, sep='\t')
        file_name = str(file)
        file_name_short = file_name[:-4]
        print(file_name_short)
        file_df.columns = ['chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2', 'score']
        file_df_input = file_df[['chrom1', 'start1', 'stop1', 'chrom2', 'start2', 'stop2']]
        file_df_input['distance'] = abs((file_df_input['start2'].astype(int)) - (file_df_input['start1'].astype(int)))
        rand_stop_list1 = get_random_contacts(file_df_input)
        rand_stop_list2 = get_random_contacts(file_df_input)
        rand_stop_list3 = get_random_contacts(file_df_input)
        file_df_input['random_stop1a'] = rand_stop_list1
        file_df_input['random_stop1a'] = file_df_input['random_stop1a'].astype(int)
        file_df_input['random_stop1b'] = file_df_input['random_stop1a'].astype(int) - file_df_input['distance'].astype(int)
        file_df_input['random_stop2a'] = rand_stop_list2
        file_df_input['random_stop2a'] = file_df_input['random_stop2a'].astype(int)
        file_df_input['random_stop2b'] = file_df_input['random_stop2a'].astype(int) - file_df_input['distance'].astype(int)
        file_df_input['random_stop3a'] = rand_stop_list3
        file_df_input['random_stop3a'] = file_df_input['random_stop3a'].astype(int)
        file_df_input['random_stop3b'] = file_df_input['random_stop3a'].astype(int) - file_df_input['distance'].astype(int)
                
        random_df1_coord1 = file_df_input[['chrom1', 'random_stop1b']]
        random_df1_coord2 = file_df_input[['chrom1', 'random_stop1a']]
        random_df1_coord1 = random_df1_coord1.reset_index(drop=True)
        random_df1_coord2 = random_df1_coord2.reset_index(drop=True)
        random_df1_coord1['random_start'] = random_df1_coord1['random_stop1b']
        random_df1_coord2['random_start'] = random_df1_coord2['random_stop1a']

        random_df2_coord1 = file_df_input[['chrom1', 'random_stop2b']]
        random_df2_coord2 = file_df_input[['chrom1', 'random_stop2a']]
        random_df2_coord1 = random_df2_coord1.reset_index(drop=True)
        random_df2_coord2 = random_df2_coord2.reset_index(drop=True)
        random_df2_coord1['random_start'] = random_df2_coord1['random_stop2b']
        random_df2_coord2['random_start'] = random_df2_coord2['random_stop2a']
        
        random_df3_coord1 = file_df_input[['chrom1', 'random_stop3b']]
        random_df3_coord2 = file_df_input[['chrom1', 'random_stop3a']]
        random_df3_coord1 = random_df3_coord1.reset_index(drop=True)
        random_df3_coord2 = random_df3_coord2.reset_index(drop=True)
        random_df3_coord1['random_start'] = random_df3_coord1['random_stop3b']
        random_df3_coord2['random_start'] = random_df3_coord2['random_stop3a']
        
        random_df1_coord1.to_csv(outpath + file_name_short + '_randomcontacts1_coord1.bed', index=None, header=None, sep='\t')
        random_df1_coord2.to_csv(outpath + file_name_short + '_randomcontacts1_coord2.bed', index=None, header=None, sep='\t')
        random_df2_coord1.to_csv(outpath + file_name_short + '_randomcontacts2_coord1.bed', index=None, header=None, sep='\t')
        random_df2_coord2.to_csv(outpath + file_name_short + '_randomcontacts2_coord2.bed', index=None, header=None, sep='\t')
        random_df3_coord1.to_csv(outpath + file_name_short + '_randomcontacts3_coord1.bed', index=None, header=None, sep='\t')
        random_df3_coord2.to_csv(outpath + file_name_short + '_randomcontacts3_coord2.bed', index=None, header=None, sep='\t')


#put randomized contacts together in one file
def get_random_aggregate_overlaps_to_file_groups(path_for_input1, path_for_input2, outpath):
    dirs1 = os.listdir(path_for_input1)
    file_list1 = []
    dirs2 = os.listdir(path_for_input2)
    file_list2 = []

    for file1 in dirs1:
        file_list1.append(file1)
    for file2 in dirs2:
        file_list2.append(file2)
    file_list1 = sorted(file_list1)
    file_list2 = sorted(file_list2)
    
    for file1, file2 in zip(file_list1, file_list2):
        file_df1 = pd.read_csv(path_for_input1 + file1, sep='\t', names=['chrom1', 'start1', 'end1'])
        file_df2 = pd.read_csv(path_for_input2 + file2, sep='\t', names=['chrom2', 'start2', 'end2'])
        file_df_together = pd.concat([file_df1, file_df2], axis=1)
        file_df_together = file_df_together.reset_index(drop=True)
        file_name = str(file1)
        file_name_short = file_name[5:-11]
        print(file_name_short)
        file_df_together.to_csv(outpath + file_name_short + '_foraggregateplot.bed', sep='\t', index=None)


