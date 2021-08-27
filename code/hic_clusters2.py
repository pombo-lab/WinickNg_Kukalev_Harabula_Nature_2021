import numpy as np, re
#import hic_tools

"""
====================================================================================================================
Python script name hic_clusters2.py
This is an additional script requried by AB-calling-v3.py
=====================================================================================================================
"""

class hicmat:
    """ Class definition for reading/analyzing HiC data sets
    """
    def __init__(self):
        self.hic_filename=''
        self.chrs=[] # name of chromosomes
        self.bins={} # coordinates of bins, currently just for the first chr.
        self.hic_data =''
        self.z_data = ''

    def read_hic_data(self, fs, skip_rows=0, skip_cols=2, order_chrs=True, chrs='', bin_col=1):
        self.hic_filename = fs
        # Extract the name of chromosomes from the filename
        if chrs=='':
            self.chrs = re.findall(r'(?<=chr).\d*',fs)
            if len(self.chrs)==1:
                self.chrs = 2*self.chrs
        else:
            self.chrs = chrs

        temp_data = np.genfromtxt(fs, skip_header=skip_rows)[:, skip_cols:]
        if order_chrs:
            with open(fs,'r') as ftemp:
                line=ftemp.readline()
                line=ftemp.readline()
                temp_chr = line.split()[0][3:]

            if temp_chr!=self.chrs[0]:
                self.chrs.reverse()
        # Set the coordinates of bins
        #self.bins[self.chrs[0]] = np.genfromtxt(fs, usecols=(bin_col), dtype=None)[skip_rows:]
        # Copy the HiC table
        self.hic_data = temp_data.copy()
        # free the memory
        del temp_data

    def clear_non_mappable(self, d1):
        N=d1.shape[0]
        mask=np.zeros(shape=d1.shape)
        # clear non-mappable data
        for i in range(0,N-1):
            for j in range(i+1, N):
                if d1[i,j]!=d1[j,i]:
                    d1[i,j] = np.nan
                    d1[j,i] = np.nan
        for i in range(0,N):
            if ((d1[i]==0).all())or(np.isnan(d1[i]).all()):
                mask[i]=True
                mask[:,i]=True
        mask[np.isnan(d1)]=True
        d1[np.isnan(d1)]=0.0
        # mask all damaged or corrupted data
        # you do not want to include those data in your analysis!
        d1=np.ma.masked_array(d1, mask=mask)
        return d1

    def calc_zscore_matrix(self, d1):
        """
        Calculate z-score matrix for an intra-chromosome Hi-C
        data set. All diagonals of the Hi-C matrix are rescaled
        to have a zero mean and standard deviation of one. So,
        all entries with the same genomic distance are rescaled
        together.

        The result has the same masked entries as the original
        Hi-C matrix.
        """
        N=d1.shape[0]
        mask=np.zeros(shape=d1.shape).astype(bool)
        ZERO = 1e-10
        # clear non-mappable data
        #for i in range(0,N-1):
        #    for j in range(i+1, N):
        #        if np.abs(d1[i,j]-d1[j,i])>ZERO:
        #            d1[i,j] = NaN
        #            d1[j,i] = NaN
        #            mask[i,j] = True
        #            mask[j,i] = True
        #for i in range(0,N):
        #    if (d1[i]==0).all():
        #        mask[i]=True
        #        mask[:,i]=True
        #mask[np.isnan(d1)]=True
        #d1[np.isnan(d1)]=0.0
        # mask all damaged or corrupted data
        # you do not want to include those data in your analysis!
        #d1=np.ma.masked_array(d1, mask=mask)
        # z-score calculation:
        #d1 = self.clear_non_mappable(d1)
        z0 = np.zeros_like(d1)
        for i in range(0,N):
            a = np.diagonal(d1, offset=i)
            rng = np.arange(N-i)
            #a = a[~a.mask]
            aprim = a[~a.mask]
            z = (a-np.mean(aprim))/np.std(aprim)
            z0[rng,rng+i] = z.copy()

        zz=z0+np.triu(z0,k=1).T
        #zz[np.isnan(zz)]=0.0
        zz = np.ma.masked_array(zz.copy(), mask=d1.mask)
        #self.hic_data = d1
        return zz

    def calc_zscore_matrix_interchromosome(self, d, nozero=False):
        d1=d.copy()
        N1=d1.shape[0]
        N2=d1.shape[1]
        mask=np.zeros_like(d1).astype(bool)
        # clear non-mappable data
        mask[np.isnan(d1)]=True
        #d1[np.isnan(d1)]=0.0
        # mask all damaged or corrupted data
        #for i in range(N1):
        #    if (d1[i]==0).all():
        #        mask[i]=True
        #for i in range(N2):
        #    if (d1[:,i]==0).all():
        #        mask[:,i]=True

        # you do not want to include those data in your analysis!
        d1=np.ma.masked_array(d1, mask=mask)
        d1 = np.ma.masked_invalid(d1)
        if nozero:
            d1=np.ma.masked_values(d1, 0)
        ########
        z0= np.ma.masked_array(np.zeros_like(d1), mask=mask)

        z0 = (d1 - d1.mean())/d1.std()

        return z0

    def calc_zscore(self, fs=''):
        if fs!='':
            self.hic_filename = fs

        if self.hic_data.size==0:
            self.read_hic_data(fs=self.hic_filename)

        if len(self.chrs)==1:
            # Intra-chromosome Hi-C
            self.z_data = self.calc_zscore_matrix(self.hic_data)
        elif len(self.chrs)==2:
            # Inter-chromosome Hi-C
            self.z_data = self.calc_zscore_matrix_interchromosome(self.hic_data)

        return self.z_data

    def save_z_matrix (self, fs=''):
        if fs=='':
            fs = self.hic_filename[:-3]+'zmat'

        N1 = self.hic_data.shape[0]
        chrs = np.array([self.chrs[0]]*N1)
        #temp_data = np.hstack((chrs.reshape(N1,1), self.bins[self.chrs[0]].reshape(N1, 1), self.z_data))
        temp_data = self.z_data
        np.savetxt(fs, temp_data, fmt='%s')
