import pandas as pd
import numpy as np
from scipy import stats
import os, sys
import glob


"""
Normalize_eigenvector_values_GAMbrain.py
Script written by: Warren Winick-Ng

This script is part of the study:
Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded by specific chromatin topologies", Nature 2021 

Date: August 16, 2021


>This script takes an input dataframe with genomic window coordinates and compartment eigenvector values.
>This script normalizes eigenvector values in each chromosome.
>Positive eigenvector values are bounded from 0 to 1
>Negative values are bounded from -1 to 0


Please refer to the Read.me document for details on usage and distribution. 

"""

def normalize_compartment_EV_per_chrom (dataframe):
    EV_list = []
    df_new = dataframe.copy()
    df_new.columns = ['coord', 'score']
    chrom_col = df_new["coord"].str.split(":", n = 1, expand = True)
    df_new['chrom'] = chrom_col[0]
    chromosomes = list (set (df_new ['chrom'].tolist()))
    for chrom in chromosomes:
        x = df_new[df_new['chrom'] == chrom] #subset chromosome
        X_sort = x.sort_values(by=['score'], ascending=False) #sort values by score
        X_A = X_sort[X_sort['score'] > 0] #subset 'A' compartment values
        X_B = X_sort[X_sort['score'] < 0] #subset 'B' compartment values
        X_A_norm = (((X_A['score'])-(X_A['score'].min()))/((X_A['score'].max())-(X_A['score'].min()))).tolist() #normalize A compartment from 0 to 1
        X_B_norm = ((X_B['score']-X_B['score'].max())/(X_B['score'].min()-X_B['score'].max())*-1).tolist() #normalize B compartment from 0 to -1
        X_final = X_A_norm + X_B_norm #create a list of A and B normalized values
        X_sort['norm_score'] = X_final
        X_sort_index = X_sort.sort_index() # put values back in the original order
        X_sort_index_list = X_sort_index['norm_score'].tolist()
        EV_list = EV_list + X_sort_index_list #add values from each chromosome to a list
    df_new['normalized_EV_values'] = EV_list
    df_new['startstop'] = chrom_col[1]
    startstop_col = df_new["startstop"].str.split("-", n = 1, expand = True)
    df_new['start'] = startstop_col[0]
    df_new['stop'] = startstop_col[1]
    df_final = df_new[['chrom', 'start', 'stop', 'normalized_EV_values']]
    return df_final

