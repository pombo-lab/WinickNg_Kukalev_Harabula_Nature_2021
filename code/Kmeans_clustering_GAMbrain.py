get_ipython().run_line_magic('matplotlib', 'inline')
import pandas as pd
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
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage


"""
Kmeans_clustering_GAMbrain.py
Script written by: Warren Winick-Ng

This script is part of the study:
Winick-Ng, W., Kukalev, A., Harabula, I., et al. "Cell-type specialization is encoded by specific chromatin topologies", Nature 2021 

Date: August 16, 2021


>This code takes an input dataframe with multiple columns for K-Means clustering and produces:
-An Elbow method graph for determining optimal cluster size
-K-means clusters fitted to the input dataframe 
-A forest dendrogram with levels based on the optimal number of clusters
-Plotting code for the resulting clusters, as a heatmap (3 example columns are shown)

This code works best if run in the Jupyter environment


Please refer to the Read.me document for details on usage and distribution. 

"""

#Visually determine the appropriate cluster number using the Elbow method
plt.figure(figsize=(10, 8))
wcss = []
for i in range(1, 11):
    kmeans = KMeans(n_clusters = i, init = 'k-means++', random_state = 42)
    kmeans.fit(input_df) #fit input data by kmeans
    wcss.append(kmeans.inertia_)
plt.plot(range(1, 11), wcss)
plt.title('The Elbow Method')
plt.xlabel('Number of clusters')
plt.ylabel('WCSS')
plt.show()


# Fitting the optimal number of K-Means clusters to the dataset
kmeans = KMeans(n_clusters = 5, init = 'k-means++', random_state = 42) # determine the optimal n_clusters from the Elbow Method plot
y_kmeans = kmeans.fit_predict(input_df)

# Begin the cluster numbering with 1 instead of 0
y_kmeans1=y_kmeans
y_kmeans1=y_kmeans+1

# New Dataframe called cluster
cluster = pd.DataFrame(y_kmeans1, index=input_df.index)
cluster.columns = ['Cluster']
#add cluster numbers to df
input_df2 = input_df.copy()
input_df2['Cluster'] = cluster['Cluster']
input_df2_sort = input_df2.sort_values('Cluster')

#get dataframe for plotting
input_df2_sort_plotting = input_df2_sort.drop(['Cluster'], axis=1)


# Hierarchical clustering for the same dataset
# Creating a dataset for hierarchical clustering
dataset2_standardized = input_df2.copy()
np.set_printoptions(precision=5, suppress=True)  # suppress scientific float notation

#Creating the linkage matrix
H_cluster = linkage(dataset2_standardized,'ward')
plt.title('Hierarchical Clustering Dendrogram (truncated)')
plt.xlabel('sample index or (cluster size)')
plt.ylabel('distance')
dendrogram(
    H_cluster,
    truncate_mode='lastp',  # show only the last p merged clusters
    p=5,  # choose the number of Kmeans clusters
    leaf_rotation=90.,
    leaf_font_size=12.,
    show_contracted=True,  # to get a distribution impression in truncated branches
)
plt.show()


#plot a figure with 3 columns clustered, min/max are set based on the values in your full dataset
plt.figure(figsize=[6,18])
plt.subplot(1,3,1)
df = input_df2_sort_plotting['name1']
sb.heatmap(df[:, np.newaxis], cmap = 'RdYlBu_r', vmin=-2.5, vmax=14.5, cbar=False, yticklabels=False, xticklabels=False)
plt.subplot(1,3,2)
df = input_df2_sort_plotting['name2']
sb.heatmap(df[:, np.newaxis], cmap = 'RdYlBu_r', vmin=-2.5, vmax=14.5, cbar=False, yticklabels=False, xticklabels=False)
plt.subplot(1,3,3)
df = input_df2_sort_plotting['name3']
sb.heatmap(df[:, np.newaxis], cmap = 'RdYlBu_r', vmin=-2.5, vmax=14.5, cbar=False, yticklabels=False, xticklabels=False)

plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)

plt.show()

