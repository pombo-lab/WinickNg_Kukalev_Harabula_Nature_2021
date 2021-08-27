#!~/.guix-profile/bin/python3

"""
====================================================================================================================
Python script name GAM.cross.contamination.test.py
Script written by Alexander Kukalev

This script is part of the study: Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded 
by specific chromatin topologies", Nature 2021

Date: August 16, 2021

This script check cross-well contaminations for the sample on the same plate by calculating Jaccard distance 
between the arrays representing the positive windows on the linear genome for every sample. It takes two input files
a) GAM QC table with a separate column indicating which samples were processed together 
b) GAM segregation table

Usage example:
>>> python3 GAM.cross.contamination.test.py $input_folder $genomic_distance 

To test the quality of genome sampling for all distances set genomic_distance to zero.
=====================================================================================================================
"""
import numpy as np
import pandas as pd
import scipy.spatial
jaccard = scipy.spatial.distance.jaccard
import matplotlib.pyplot as plt
import seaborn as sns

def convert_to_jaccard_dists(seg_table, seg_table_cols):    
    values = []
    for i in range(len(list(seg_table.columns))):
        for j in range(len(list(seg_table.columns))):
            values.append(jaccard(np.array(seg_table.iloc[:,i]),np.array(seg_table.iloc[:,j])))
        
    values = np.array(values).reshape(len(seg_table.columns),len(seg_table.columns))
    values_df = (pd.DataFrame(values))
    values_df.columns = seg_table_cols
    values_df['index'] = seg_table_cols
    values_df.set_index('index',inplace=True)
    return values_df

def process_results(values_df,seg_table_cols,similarity_threshold):
    fig = sns.clustermap(values_df)
    new_order = []
    for i in fig.dendrogram_row.reordered_ind:
        new_order.append(seg_table_cols[i])
    clustered_values_df = values_df[new_order]
    clustered_values_df = clustered_values_df.reindex(new_order)
    boolean_df = (clustered_values_df > similarity_threshold).astype(int)
    for i in range(len(boolean_df)):
        boolean_df.iloc[i,i] = 1
    bad_wells = []
    for i in seg_table_cols:
        if 0 in list(boolean_df.loc[i]):
            bad_wells.append(i)
    plt.close('all')
    return (bad_wells)

def contamination_bias_wrapper(subset_seg_table, similarity_threshold):
    seg_table_cols = subset_seg_table.columns
    values_df = convert_to_jaccard_dists(subset_seg_table, seg_table_cols)
    cross_well_contaminated_samples = process_results(values_df, seg_table_cols, similarity_threshold)
    return cross_well_contaminated_samples

def cross_well_contamination_per_plate (input_stats_df, segregation_df, plate_column, similarity_threshold = 0.4):
    # This is the main function
    stats_df = input_stats_df.copy()
    list_of_plates = set (stats_df[plate_column].tolist())
    for plate in list_of_plates:
        subset_plate = stats_df.query (f'{plate_column}==@plate')
        samples_on_plate = subset_plate ['Simple_name'].tolist()
        segregation_per_plate = segregation_df [samples_on_plate]
        cross_well_contaminated_samples = contamination_bias_wrapper(segregation_per_plate, similarity_threshold)
        if cross_well_contaminated_samples != []:
            for sample in cross_well_contaminated_samples:
                sample_QC_state = stats_df.query ('Simple_name == @sample')['QC'].tolist()
                if sample_QC_state == ['passed']:
                    print ('Plate --> ' + plate)
                    print ('Sample is cross-contaminated ' + sample)
                    stats_df.loc [stats_df['Simple_name'] == sample, 'QC'] = 'not_passed'
                    stats_df.loc [stats_df['Simple_name'] == sample, 'QC_reason'] = 'cross_contaminated'
    return stats_df

seg_table = pd.read_csv ('path/to/GAM/segregation/table', sep='\t', index_col=[0,1,2])
GAM_stats_df = pd.read_csv ('path/to/GAM/QC/stats/table', sep='\t')
GAM_QC_samples = cross_well_contamination_per_plate (GAM_stats_df, seg_table, 'Plate')