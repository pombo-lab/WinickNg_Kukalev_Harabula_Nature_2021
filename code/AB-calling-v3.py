#! /usr/bin/env python3

"""
====================================================================================================================
Python script name AB-calling-v3.py
Script written by Ibai Irastorza-Azcarate, original version written by Ehsan Irani 

This script is part of the study: Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded 
by specific chromatin topologies", Nature 2021

Date: August 16, 2021

This script performs compartment calling on GAM data

Usage example:
>>> python3 AB-calling-v3.py -i $input_cosegregation_matrix --output_folder folder_name --output_file file_name --corr-bedfile GC_content_bed --chrom chrom_name --clear-nonmap True

=====================================================================================================================
"""

import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA as sklearnPCA
#from sklearn import PCA as sklearnPCA
import hic_clusters2 as hic
import argparse
import pybedtools
from scipy.stats import pearsonr
import os
# %%
######################################################
def divide_over_expected (d1):
    N=d1.shape[0]
    mask=np.zeros(shape=d1.shape).astype(bool)
    ZERO = 1e-10
    z0 = np.zeros_like(d1)
    for i in range(0,N):
        a = np.diagonal(d1, offset=i)
        rng = np.arange(N-i)
        #a = a[~a.mask]
        aprim = a[~a.mask]
        aprim[np.isnan(aprim)] = 1
        aprim_mean = np.mean(aprim)
        z = a/aprim_mean
        z0[rng,rng+i] = z.copy()

    zz = z0+np.triu(z0,k=1).T
    zz = np.ma.masked_array(zz.copy(), mask=d1.mask)
    return zz
######################################################
# %%

parser = argparse.ArgumentParser()
group_gc = parser.add_argument_group()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--input', help='input tabular file', type=str, dest='input', required=True)
requiredNamed.add_argument('--output_folder', help='output directory', type=str, dest='folder_name', required=True)
parser.add_argument('--output_file', help='output bed file', type=str, dest='output', default='AB-compartment.bed')
group_gc.add_argument('--corr-bedfile',\
help='bedfile to correlate and indicate the A compartment', type=str, dest='corr_bedfile', default=None)
group_gc.add_argument('--chrom', help='chrom in the corr-bedfile', default=None, type=str, dest='chrom')
parser.add_argument('--clear-nonmap',\
help='pass true if you want to clear and mask non-mappable areas int the input file', type=bool, default=False, dest='clear_non_map')
#pareser.add_argument('--clear-nonmap', type=bool, default=True, dest='clear_nonmap')
parser.add_argument('--skip-rows', '--sr', help='number of rows to skip in reading the input file', type=int, default=0, dest='skip_rows')
parser.add_argument('--skip-cols', '--sc', help='number of columns to skip in reading the input file', type=int, default=0, dest='skip_cols')
parser.add_argument('--n-comp', '--nc', help='number of components to write in the output file', type=int, default=1, dest='nc')
parser.add_argument('--style',  help='1 (observed/expected), 2 (zscore), 3 (zscore of log)', type=int, default=1, dest='style')
args = parser.parse_args()

if not os.path.exists(args.folder_name):
    os.mkdir(args.folder_name)
    print("Directory ", args.folder_name, " created")
else:
    print("Directory ", args.folder_name, " already exists")

# %%
########## Read the inout matrix ##########
hic1 = hic.hicmat()
hic1.read_hic_data(args.input, skip_rows=args.skip_rows, skip_cols=args.skip_cols, chrs=args.chrom, order_chrs=False)
if args.clear_non_map:
    hic1.hic_data = hic1.clear_non_mappable(hic1.hic_data)
########## Read the inout matrix ##########
# Create the expected matrix
# create the observed over expected
K = pd.DataFrame(divide_over_expected(hic1.hic_data))
#K.to_csv("test.txt",sep="\t")


#keep lines/rows that contain at least 90% of non-NANs
#thresh = len(K.index)*0.9
#K = K.dropna(axis='rows',thresh = thresh) #drop NANs in both rows and cols.
#remaining_rows = K.index
#K = K.dropna(axis='columns',thresh = thresh) #drop NANs in both rows and cols.


K = K.mask(np.isnan(K.values), 0)
data_unfiltered = K
# count zeros
zero_count = (K==0).sum(axis=1)
remaining_rows = zero_count[zero_count < (len(K.index)*0.25)] # at least 75% of non-zeros
remaining_rows = remaining_rows.index.values
K = K[remaining_rows]
K = K.iloc[remaining_rows]

#K = divide_over_expected(hic1.hic_data)
# create the correlation save_z_matrix

if args.style==1:
    data = K
elif args.style==2:
    Z = hic1.calc_zscore_matrix(hic1.hic_data)
    data = pd.DataFrame(Z.data)
    data[Z.mask] = np.NaN
elif args.style==3:
    log_mat = hic1.hic_data+1
    log_Z = hic1.calc_zscore_matrix(np.log10(log_mat))
    data = pd.DataFrame(log_Z)
    data[log_Z.mask] = np.NaN

print('Calculate the correlation matrix of Observed/Expected...\n')
corr_mat = data.corr()
corr_mat = corr_mat.mask(np.isnan(corr_mat.values), 0)
data_unfiltered = data_unfiltered.corr()
data_unfiltered = data_unfiltered.mask(np.isnan(data_unfiltered.values),0)

# %%
###########################################
sklearn_pca = sklearnPCA(n_components = args.nc, svd_solver='full')
corr_sign = 1


# DO WE WANT TO CORRELATE DATA OR DATA2???
y_sklearn = sklearn_pca.fit_transform(corr_mat)

if args.corr_bedfile!=None:
    corr_bedfile = pybedtools.BedTool(args.corr_bedfile)
    df_corr_bedfile = corr_bedfile.to_dataframe()
    if args.chrom!=None:
        df_corr_bedfile = df_corr_bedfile[df_corr_bedfile['chrom']==args.chrom]

    corr_data = df_corr_bedfile[df_corr_bedfile.columns[-1]]
df_corr_bedfile = df_corr_bedfile.iloc[remaining_rows]
corr_data = corr_data.iloc[remaining_rows]
# %%
# ### write the output
out_columns = ['chrom', 'start', 'end', 'AB']
out_columns = out_columns + ['score%d' %i for i in range(args.nc)]
if args.chrom==None:
    out_columns = out_columns[1:]

output_df = pd.DataFrame(columns=out_columns)
mask = np.diagonal(hic1.hic_data.mask)
mask = [x for n,x in enumerate(mask) if n in remaining_rows]
mask = np.array(mask)
output_df.start = df_corr_bedfile.start[~mask]
output_df.end = df_corr_bedfile.end[~mask]
header_str = '# n_pca_comp %di\n'%args.nc
# ### Other PCAs
pca_id = 0

""" Change to 0 to 3 for all PCA """

#for i in range(0, args.nc): # IF we want to get the PCA depending on GC, uncomment this
for i in range(0, 3): #now we always take the first component as Compartments
    temp_corr = pearsonr(corr_data, y_sklearn[:,i])[0]
    print ('PCA {} comp. corr.: {}'.format(i, temp_corr)) 
    header_str = header_str + '#PCA {} comp. corr.: {}\n'.format(i, temp_corr)
    if i==0:
        temp_corr0 = temp_corr
    elif np.abs(temp_corr)>np.abs(temp_corr0):
        pca_id = i
        temp_corr0 = temp_corr

    corr_sign = np.sign(temp_corr)
    output_df['score%d'%i] = corr_sign*y_sklearn[:,i][~mask]
# ### First PCA
#y_sklearn = corr_sign * y_sklearn
corr_sign = np.sign(temp_corr0)
print (pca_id, temp_corr0, corr_sign)
A_id = corr_sign*y_sklearn[:,pca_id]>=0
B_id = corr_sign*y_sklearn[:,pca_id]<0
AB = np.array(['A']*len(mask))
AB[B_id] = 'B'
output_df.AB = AB[~mask]
if args.chrom!=None:
    output_df.chrom = args.chrom
output_df.to_csv(args.folder_name + "/" + args.output, sep='\t', header=header_str, index=False)
style_name = {
1 : "obs_over_exp",
2 : "zscore",
3 : "yscore_log"
}
corr_mat.to_csv("{}/{}_corr_{}".format(args.folder_name,style_name[args.style], args.chrom), sep='\t', header=header_str, index=False)
data.to_csv("{}/{}_input_{}".format(args.folder_name,style_name[args.style], args.chrom), sep='\t', header=header_str, index=False)
data_unfiltered.to_csv("{}/{}_corr_unfiltered_{}".format(args.folder_name,style_name[args.style], args.chrom), sep='\t', header=header_str, index=False)
