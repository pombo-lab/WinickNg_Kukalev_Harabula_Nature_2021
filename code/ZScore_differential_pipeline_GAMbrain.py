import sys
import glob
import gzip
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import norm
import pysam
import numpy as np

"""
ZScore_differential_pipeline_GAMbrain.py
Script written by: Warren Winick-Ng

This script is part of the study:
Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded by specific chromatin topologies", Nature 2021 

Date: August 16, 2021


This script takes two input GAM segregation tables, and will generate for each chromosome:
-NPMI matrices 
-Z-Score matrices
-Differential Z-Score matrices
-Top X% significant differential matrices
For each of the above, it is possible to produce full length and/or distance cut-off matrices (cutoff parameter)

This script also determines lowly detected windows (default 2%) in each input segregation table and removes them from both datasets


The main function is 'Segregation_to_ZScore_differential';
:param ID_1,ID_2: dataset ID for naming the output file
:param seg_table1: path to first input matrix or open python file
:param seg_table2: path to second input matrix or open python file
:param region: chromosome in format 'chr#'
:param cutoff: provide a distance cutoff for output matrices in Mb; input 0 to skip cutoff
:param threshold: provide the threshold percentage for top X% differential contacts
:param outfile: path to output directory

Usage example:
>>> ./ZScore_differential_pipeline_GAMbrain.py 'ID_1' 'ID_2' path_to_seg_table1 path_to_seg_table2 'chr1' 5000000 5 path_to_output_folder

Please refer to the Read.me document for details on usage and distribution. 

"""

def open_segregation(path_or_buffer):
	"""
	Open a segregation table from a file.
	:param path_or_buffer: Path to imput segregation table, or open python file object.
	:returns: :ref:`segregation table <segregation_table>`
	"""
	open_seg = pd.read_csv(path_or_buffer, sep='\t') #, nrows=100
	open_seg.set_index(['chrom','start','stop'], inplace = True)	
	return open_seg

def calculate_NPMI(seg_matrix):
	"""
	Calculate NPMI matrix from a segregation table
	"""
	
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

def windows_detection_frequencies (segregation_table):
	"""
	Calculate window detection frequencies for each window in a chromosome from a segregation table
	"""
	if len (segregation_table.index.names) == 1: # Check if supplied segregation table has multiindex
		segregation_table.set_index(['chrom','start','stop'], inplace = True)
	number_of_samples = len(list(segregation_table))
	frequencies_list = []
	for index, row in segregation_table.iterrows(): # Itterate over rows of the segregation table
		number_of_positive = (row == 1).values.sum()
		frequency = float(number_of_positive/number_of_samples)
		frequencies_list.append(frequency)
	frequencies_list = np.around(np.array(frequencies_list),4) # Keep first 4 digits after decimal point
	frequencies_df = pd.DataFrame(index=segregation_table.index)
	frequencies_df['data'] = pd.DataFrame(frequencies_list, index=segregation_table.index)

	return frequencies_df

def get_bad_index_list_from_WDF(frequencies_df1, frequencies_df2):
	"""
	Generate a list of windows which need to be masked in the Z-score calculation.
	:param frequencies_df1, frequencies_df2: dataframe of window detection frequencies with column 'data' 
	"""
	
	#generate df for each frequency list that has zeros in place of lowest 2% of frequencies
	noZerodf1 = frequencies_df1.loc[~(frequencies_df1==0).all(axis=1)]
	greater_than_2per_array1 = np.where(frequencies_df1['data'] < (np.nanpercentile(noZerodf1, 2)), 0, frequencies_df1['data']) #determine lowest 2% of windows from dataframe1
	greater_than_2per_array_df1 = pd.DataFrame(greater_than_2per_array1, index=frequencies_df1.index)
	greater_than_2per_array_df1.columns = ['data']
	noZerodf2 = frequencies_df2.loc[~(frequencies_df2==0).all(axis=1)]
	greater_than_2per_array2 = np.where(frequencies_df2['data'] < (np.nanpercentile(noZerodf2, 2)), 0, frequencies_df2['data']) #determine lowest 2% of windows from dataframe2
	greater_than_2per_array_df2 = pd.DataFrame(greater_than_2per_array2, index=frequencies_df2.index)
	greater_than_2per_array_df2.columns = ['data']
	
	#generate list of zero windows for each dataset
	bad_index_zeros1 = greater_than_2per_array_df1.loc[greater_than_2per_array_df1['data'] == 0]
	bad_index_zeros2 = greater_than_2per_array_df2.loc[greater_than_2per_array_df2['data'] == 0]
	bad_index_zeros1_reset = bad_index_zeros1.reset_index()
	bad_index_zeros2_reset = bad_index_zeros2.reset_index()
	bad_index_zeros1_reset['index_list'] = (bad_index_zeros1_reset['chrom'] + ':' + bad_index_zeros1_reset['start'].map(str) + '-' + bad_index_zeros1_reset['stop'].map(str)).astype(str)
	bad_index_zeros2_reset['index_list'] = (bad_index_zeros2_reset['chrom'] + ':' + bad_index_zeros2_reset['start'].map(str) + '-' + bad_index_zeros2_reset['stop'].map(str)).astype(str)
	zeros_index_list1 = bad_index_zeros1_reset['index_list'].values.tolist()
	zeros_index_list2 = bad_index_zeros2_reset['index_list'].values.tolist()
	
	#determine non-overlapping <2% windows in both datasets and create a list with only with these windows
	greater_than_2per_array_df1['match'] = np.where((greater_than_2per_array_df2['data'] == 0) == (greater_than_2per_array_df1['data'] > 0), 'Yes', 'No') 
	greater_than_2per_array_df2['match'] = np.where((greater_than_2per_array_df1['data'] == 0) == (greater_than_2per_array_df2['data'] > 0), 'Yes', 'No') 
	bad_sample_df1_noNA = (greater_than_2per_array_df1.where(greater_than_2per_array_df1['match'].isin(['Yes']))).dropna()
	bad_sample_df2_noNA = (greater_than_2per_array_df2.where(greater_than_2per_array_df2['match'].isin(['Yes']))).dropna()
	bad_sample_df1_noNA_reset = bad_sample_df1_noNA.reset_index()
	bad_sample_df2_noNA_reset = bad_sample_df2_noNA.reset_index()
	bad_sample_df1_noNA_reset['index_list'] = (bad_sample_df1_noNA_reset['chrom'] + ':' + bad_sample_df1_noNA_reset['start'].map(str) + '-' + bad_sample_df1_noNA_reset['stop'].map(str)).astype(str)
	bad_sample_df2_noNA_reset['index_list'] = (bad_sample_df2_noNA_reset['chrom'] + ':' + bad_sample_df2_noNA_reset['start'].map(str) + '-' + bad_sample_df2_noNA_reset['stop'].map(str)).astype(str)
	bad_index_list1 = bad_sample_df1_noNA_reset['index_list'].values.tolist()
	bad_index_list2 = bad_sample_df2_noNA_reset['index_list'].values.tolist()
	
	#create a list of windows to be masked for each df
	final_indexlist1 = bad_index_list1 + zeros_index_list1
	final_indexlist2 = bad_index_list2 + zeros_index_list2
	
	return final_indexlist1, final_indexlist2
	
	
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
	return Z_matrix


def get_top5_bottom5_matrices (diff_matrix, threshold):
	"""
	Calculate top X% of contacts for each dataset in a differential Z-Score matrix
	"""
	flat_matrix = np.matrix.flatten(diff_matrix)
	data_flat_matrix = flat_matrix[~np.isnan(flat_matrix)]
	threshold_int = int(threshold)
	threshold_1 = threshold_int/100
	threshold_2 = 1 - threshold_1
	print('theshold1 = ' + str(threshold_1))
	print('theshold2 = ' + str(threshold_2))
	(mu_diff_matrix, sigma_diff_matrix) = norm.fit(data_flat_matrix)
	top5_matrix2 = norm.ppf(threshold_1, loc=mu_diff_matrix, scale=sigma_diff_matrix)
	top5_matrix1 = norm.ppf(threshold_2, loc=mu_diff_matrix, scale=sigma_diff_matrix)
	matrix1 = diff_matrix.copy()
	matrix2 = diff_matrix.copy()
	matrix1[matrix1 <= top5_matrix1] = np.nan
	matrix2[matrix2 > top5_matrix2] = np.nan
	return matrix1, matrix2


def get_cutoff_matrix(input_matrix, chrom, cutoff, segregation_chrom):
	"""
	Restrict the distance of the matrix to a specified cutoff
	"""
	chrom_sizes = {"chr1" : 195471971,"chr2" : 182113224,"chr3" : 160039680,"chr4" : 156508116, \
		"chr5" : 151834684,"chr6" : 149736546,"chr7" : 145441459,"chr8" : 129401213, \
		"chr9" : 124595110,"chr10" : 130694993, "chr11" : 122082543,"chr12" : 120129022, \
		"chr13" : 120421639,"chr14" : 124902244,"chr15" : 104043685,"chr16" : 98207768, \
		"chr17" : 94987271,"chr18" : 90702639,"chr19" : 61431566,"chrX": 171031299,"chrY" : 91744698}  
	value_list = []
	end_region = int(chrom_sizes[chrom])
	col1_name = segregation_chrom.columns[0]
	col2_name = segregation_chrom.columns[1]
	coordinate1 = col1_name.split(':')[1]
	coordinate2 = col2_name.split(':')[1]
	start1 = int(coordinate1.split('-')[0])
	start2 = int(coordinate2.split('-')[0])
	interval_size = int(start2 - start1)
	cut_off = int(int(cutoff)/int(interval_size))
	matrix_10Mb = np.ma.masked_invalid(input_matrix)
	for i in range(cut_off+1, len(matrix_10Mb)):
		diag = np.diagonal(matrix_10Mb, offset=i)
		fillNA_diag = np.nan 
		np.fill_diagonal(matrix_10Mb[i:,:], fillNA_diag)
		np.fill_diagonal(matrix_10Mb[:,i:], fillNA_diag)
	return matrix_10Mb

  
def Segregation_to_ZScore_differential(argv):
	"""
	Generate differential Z-Score matrices for a whole chromosome.
	:param ID_1,ID_2: dataset ID for naming the output file
	:param seg_table1: path to first input matrix or open python file
	:param seg_table2: path to second input matrix or open python file
	:param region: chromosome in format 'chr#'
	:param cutoff: provide a distance cutoff for output matrices in Mb, input 0 to skip
	:param threshold: provide the threshold percentage for top X% differential contacts
	:param outfile: path to output directory
	"""	
	if len(argv)==8:
		(ID_1, ID_2, seg_table1, seg_table2, region, cutoff, threshold, outfile) = argv
	else:
		print('Wrong number of arguments')
		sys.exit()
	
	#open segregation tables and subset chromosome
	open1 = open_segregation(seg_table1)
	open2 = open_segregation(seg_table2)
	segregation_chrom1 = open1.loc[region]
	segregation_chrom2 = open2.loc[region]
	segregation_matrix1 = segregation_chrom1.values
	segregation_matrix2 = segregation_chrom2.values
	
	#create bin names [chr:start-end]
	cnames = segregation_chrom1.index.to_frame(index=None)
	cnames["bin"]=region
	bins=cnames[("bin")].str.cat(cnames[("start")].astype(str), sep=":").str.cat(cnames[("stop")].astype(str), sep="-")

	#generate NPMI matrix for chromosome
	NPMI_matrix1 = calculate_NPMI(segregation_matrix1)
	NPMI_matrix2 = calculate_NPMI(segregation_matrix2)
	
	#generate list of windows to be masked in Z-Score matrix
	WDF_df1 = windows_detection_frequencies(open1)
	WDF_df2 = windows_detection_frequencies(open2)
	mask_indexlist1, mask_indexlist2 = get_bad_index_list_from_WDF(WDF_df1, WDF_df2)
	mask_indexlist1_chrom = [k for k in mask_indexlist1 if region+':' in k]
	mask_indexlist2_chrom = [k for k in mask_indexlist2 if region+':' in k]
	
	#set indices to be masked in Z-score matrices to NaN
	NPMI_dataframe1_formasking = pd.DataFrame(NPMI_matrix1)
	NPMI_dataframe2_formasking = pd.DataFrame(NPMI_matrix2)
	NPMI_dataframe1_formasking.columns=bins.astype(str)
	NPMI_dataframe2_formasking.columns=bins.astype(str)
	NPMI_dataframe1_formasking_final = NPMI_dataframe1_formasking.rename(index=bins.astype(str))
	NPMI_dataframe2_formasking_final = NPMI_dataframe2_formasking.rename(index=bins.astype(str))
	NPMI_dataframe1_formasking_final.loc[mask_indexlist1_chrom, :] = np.nan
	NPMI_dataframe2_formasking_final.loc[mask_indexlist2_chrom, :] = np.nan
	NPMI_dataframe1_formasking_final.loc[:, mask_indexlist1_chrom] = np.nan
	NPMI_dataframe2_formasking_final.loc[:, mask_indexlist2_chrom] = np.nan
	NPMI_masked1 = NPMI_dataframe1_formasking_final.values
	NPMI_masked2 = NPMI_dataframe2_formasking_final.values
		
	#generate ZScores from NPMI and take differential
	Zmatrix1 = zscore_matrix(NPMI_masked1)
	Zmatrix2 = zscore_matrix(NPMI_masked2)
	z_diff = Zmatrix1 - Zmatrix2

	#generate top X% differential Z-Scores for each dataset
	top5_z_diff_Zmatrix1, top5_z_diff_Zmatrix2 = get_top5_bottom5_matrices(z_diff, threshold)

	#create dataframes for all of generated matrices
	NPMI_dataframe1 = pd.DataFrame(NPMI_matrix1)
	NPMI_dataframe2 = pd.DataFrame(NPMI_matrix2)
	Zmatrix1_dataframe = pd.DataFrame(Zmatrix1)
	Zmatrix2_dataframe = pd.DataFrame(Zmatrix2)
	z_diff_dataframe = pd.DataFrame(z_diff)
	top5_z_diff_Zmatrix1_dataframe = pd.DataFrame(top5_z_diff_Zmatrix1)
	top5_z_diff_Zmatrix2_dataframe = pd.DataFrame(top5_z_diff_Zmatrix2)

	#rename columns/index NPMI, ZScore, differential, top matrix files	
	NPMI_dataframe1_formasking_final.columns=bins.astype(str)
	NPMI_dataframe2_formasking_final.columns=bins.astype(str)
	Zmatrix1_dataframe.columns=bins.astype(str)
	Zmatrix2_dataframe.columns=bins.astype(str)	
	z_diff_dataframe.columns=bins.astype(str)
	top5_z_diff_Zmatrix1_dataframe.columns=bins.astype(str)
	top5_z_diff_Zmatrix2_dataframe.columns=bins.astype(str)
	NPMI_dataframe1_final = NPMI_dataframe1_formasking_final.rename(index=bins.astype(str))
	NPMI_dataframe2_final = NPMI_dataframe2_formasking_final.rename(index=bins.astype(str))
	Zmatrix1_final = Zmatrix1_dataframe.rename(index=bins.astype(str))
	Zmatrix2_final = Zmatrix2_dataframe.rename(index=bins.astype(str))
	z_diff_final = z_diff_dataframe.rename(index=bins.astype(str))
	top5_z_diff_Zmatrix1_dataframe_final = top5_z_diff_Zmatrix1_dataframe.rename(index=bins.astype(str))
	top5_z_diff_Zmatrix2_dataframe_final = top5_z_diff_Zmatrix2_dataframe.rename(index=bins.astype(str))

	#save all desired files
	NPMI_dataframe1_final.to_csv(outfile + str(ID_1) + '_' + region + '_NPMI.txt.gz',
		sep="\t",
		index=True,
		header=True,
		compression='gzip')	
	NPMI_dataframe2_final.to_csv(outfile + str(ID_2) + '_' + region + '_NPMI.txt.gz',
		sep="\t",
		index=True,
		header=True,
		compression='gzip')	
	Zmatrix1_final.to_csv(outfile + str(ID_1) + '_' + region + '_ZScore.txt.gz',
		sep="\t",
		index=True,
		header=True,
		compression='gzip')	
	Zmatrix2_final.to_csv(outfile + str(ID_2) + '_' + region + '_ZScore.txt.gz',
		sep="\t",
		index=True,
		header=True,
		compression='gzip')
	z_diff_final.to_csv(outfile + str(ID_1) + '_' + str(ID_2) + '_' + region + '_ZScore_diff.txt.gz',
		sep="\t",
		index=True,
		header=True,
		compression='gzip')
	top5_z_diff_Zmatrix1_dataframe_final.to_csv(outfile + str(ID_1) + '_' + str(ID_2) + '_' + region + '_ZScore_diff_top' + str(threshold) + '%_' + str(ID_1) + '.txt.gz',
		sep="\t",
		index=True,
		header=True,
		compression='gzip')
	top5_z_diff_Zmatrix2_dataframe_final.to_csv(outfile + str(ID_1) + '_' + str(ID_2) + '_' + region + '_ZScore_diff_top' + str(threshold) + '%_' + str(ID_2) + '.txt.gz',
		sep="\t",
		index=True,
		header=True,
		compression='gzip')

	#Mask regions above input cutoff distance and create dataframes
	if int(cutoff) > 0:
		z_diff_final_short = get_cutoff_matrix(z_diff, region, cutoff, NPMI_dataframe1_final)
		top5_z_diff_Zmatrix1_final_short = get_cutoff_matrix(top5_z_diff_Zmatrix1, region, cutoff, NPMI_dataframe1_final)
		top5_z_diff_Zmatrix2_final_short = get_cutoff_matrix(top5_z_diff_Zmatrix2, region, cutoff, NPMI_dataframe1_final)
		NPMImatrix1_short = get_cutoff_matrix(NPMI_matrix1, region, cutoff, NPMI_dataframe1_final)
		NPMImatrix2_short = get_cutoff_matrix(NPMI_matrix2, region, cutoff, NPMI_dataframe1_final)
		Zmatrix1_short = get_cutoff_matrix(Zmatrix1, region, cutoff, NPMI_dataframe1_final)
		Zmatrix2_short = get_cutoff_matrix(Zmatrix2, region, cutoff, NPMI_dataframe1_final)
		z_diff_short_dataframe = pd.DataFrame(z_diff_final_short)
		top5_z_diff_Zmatrix1_short_dataframe = pd.DataFrame(top5_z_diff_Zmatrix1_final_short)
		top5_z_diff_Zmatrix2_short_dataframe = pd.DataFrame(top5_z_diff_Zmatrix2_final_short)
		NPMImatrix1_short_dataframe = pd.DataFrame(NPMImatrix1_short)
		NPMImatrix2_short_dataframe = pd.DataFrame(NPMImatrix2_short)
		Zmatrix1_short_dataframe = pd.DataFrame(Zmatrix1_short)
		Zmatrix2_short_dataframe = pd.DataFrame(Zmatrix2_short)
		z_diff_short_dataframe.columns=bins.astype(str)
		
		#rename columns/index NPMI, ZScore, differential, top cutoff matrix files
		top5_z_diff_Zmatrix1_short_dataframe.columns=bins.astype(str)
		top5_z_diff_Zmatrix2_short_dataframe.columns=bins.astype(str)
		NPMImatrix1_short_dataframe.columns=bins.astype(str)
		NPMImatrix2_short_dataframe.columns=bins.astype(str)
		Zmatrix1_short_dataframe.columns=bins.astype(str)
		Zmatrix2_short_dataframe.columns=bins.astype(str)
		z_diff_short_dataframe_final = z_diff_short_dataframe.rename(index=bins.astype(str))
		top5_z_diff_Zmatrix1_short_dataframe_final = top5_z_diff_Zmatrix1_short_dataframe.rename(index=bins.astype(str))
		top5_z_diff_Zmatrix2_short_dataframe_final = top5_z_diff_Zmatrix2_short_dataframe.rename(index=bins.astype(str))
		NPMImatrix1_short_dataframe_final = NPMImatrix1_short_dataframe.rename(index=bins.astype(str))
		NPMImatrix2_short_dataframe_final = NPMImatrix2_short_dataframe.rename(index=bins.astype(str))
		Zmatrix1_short_dataframe_final = Zmatrix1_short_dataframe.rename(index=bins.astype(str))
		Zmatrix2_short_dataframe_final = Zmatrix2_short_dataframe.rename(index=bins.astype(str))
		
		#save cutoff matrices
		z_diff_short_dataframe_final.to_csv(outfile + str(ID_1) + '_' + str(ID_2) + '_' + region + '_ZScore_diff_' + cutoff + '_cutoff.txt.gz', sep="\t", index=True, header=True, compression='gzip')
		top5_z_diff_Zmatrix1_short_dataframe_final.to_csv(outfile + str(ID_1) + '_' + str(ID_2) + '_' + region + '_ZScore_diff_top' + str(threshold) + 'per_' + cutoff + '_cutoff_' + str(ID_1) + '.txt.gz', sep="\t", index=True, header=True, compression='gzip')
		top5_z_diff_Zmatrix2_short_dataframe_final.to_csv(outfile + str(ID_1) + '_' + str(ID_2) + '_' + region + '_ZScore_diff_top' + str(threshold) + 'per_'+  cutoff + '_cutoff_' + str(ID_2) + '.txt.gz', sep="\t", index=True, header=True,	compression='gzip')
		NPMImatrix1_short_dataframe_final.to_csv(outfile + str(ID_1) + '_' + region + '_NPMI_' + cutoff + '_cutoff.txt.gz', sep="\t", index=True, header=True,	compression='gzip')
		NPMImatrix2_short_dataframe_final.to_csv(outfile + str(ID_2) + '_' + region + '_NPMI_' + cutoff + '_cutoff.txt.gz', sep="\t", index=True, header=True, compression='gzip')
		Zmatrix1_short_dataframe_final.to_csv(outfile + str(ID_1) + '_' + region + '_ZScore_' + cutoff + '_cutoff.txt.gz', sep="\t", index=True, header=True, compression='gzip')
		Zmatrix2_short_dataframe_final.to_csv(outfile + str(ID_2) + '_' + region + '_ZScore_' + cutoff + '_cutoff.txt.gz',	sep="\t", index=True, header=True, compression='gzip')
	
	elif int(cutoff) == 0:
		print('Skipping cutoff matrices for ' + str(region))
		pass

	print(region + ' done!')
	pass

#%%

Segregation_to_ZScore_differential(sys.argv[1:])
