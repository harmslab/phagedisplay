import pyximport
pyximport.install()
import DistMatrix

from sklearn.cluster import *
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform

import numpy as np
import matplotlib.pyplot as plt


class Cluster():
    """
    Performs a given clustering algorithm on a distance matrix.
    args
        dist_matrix: distance matrix to perform clustering on
        cluster_alg: clustering algorithm to use on matrix
        factor: for setting a parameter on calculating clusters/number of predicted clusters. (epsilon DBSCAN)
        num: DBSCAN - min_samples
    """
    
    def __init__(self, dist_matrix, cluster_alg, factor = None, num = None):
        self._dist_matrix = dist_matrix
        self._cluster_alg = cluster_alg
        self._factor = factor
        self._num = num

    def freqHist(self):
        """
        make histogram plot of frequency of sequence distance scores in the given distance matrix.
        """
        
        plt.figure();
        self._dist_matrix.plot(kind='hist', legend=False, orientation='horizontal')

    def cluster(self):
        """
        cluster using DBSCAN or hierarchal clustering. DBSCAN works best on larger datasets
        while hierarchal works best on smaller datasets.
        DBSCAN - 
            min_samples: the minimum number of samples for the point to considered to be a central cluster point.
                a higher number will give a smaller number of large clusters with a large noise (-1) cluster while
                a lower number will give a large number of small clusters with a very small or no noise cluster.
            
            epsilon: use of freqHist can help determine a factor to use, max distance radius around a point to grab
                cluster points from. 
        """

        X = self._dist_matrix

        # compute scipy hierarchal clustering
        if self._cluster_alg == 'hierarchal':
            condensed = squareform(X)
            linkage = hcl.average(condensed)
            clusters = hcl.fcluster(linkage, self._factor, criterion = 'distance')

        # compute DBSCAN clustering
        elif self._cluster_alg == 'DBSCAN':
            db = DBSCAN(eps = self._factor, min_samples = self._num, metric = 'precomputed', 
                       leaf_size = 300, algorithm = 'ball_tree').fit(X.as_matrix() )
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            clusters = db.labels_

        else:
            raise ValueError('algorithm unavailable')

        # get number of clusters
        n_clusters = len(np.unique(clusters)) - (1 if -1 in clusters else 0)
        cluster_labels = pd.DataFrame({'Sequences' : X.index, 
                                       'Cluster' : clusters})
        
    
        return cluster_labels


class Membership():
    """
    accessing cluster membership of given cluster in a pandas dataframe format.
    """
    
    def __init__(self, cluster):
        self._cluster = cluster.groupby('Cluster')
        
    def getAll(self):
        """
        return dictionary of all clusters and member IDs. 
        """
        
        return self._cluster.groups
    
    def getCluster(self, clust_num):
        """
        return specified cluster and members.
        """
        
        return self._cluster.get_group(clust_num)
    
    def getCount(self):
        """
        return count of size of each cluster.
        """
        
        return self._cluster.count()
    
    def getOverlap(self, x, y, clust2):
        """
        compare two clusters to see amount of overlap.
        """
        
        overlap = []

        for i in range(x):
            for j in range(y):
                a = self.getCluster(i)
                b = clust2.getCluster(j)
                percent_overlap = round(((len(np.intersect1d(a.as_matrix(['Sequences']), b.as_matrix(['Sequences'])))
                                         /len(a))*100), 2)

                # returns non-unique values in both arrays.
                overlap.append((i, j, len(a), len(b), percent_overlap))

        col = ['cluster a', 'cluster b', 'length a', 'length b', '% overlap']
        pd.DataFrame(overlap, columns = col)
        
        return overlap
    
    def toCSV(self, cluster, file):
        """
        save a cluster to a csv file.
        """
        
        self.getCluster(cluster)['Sequences'].to_csv(file, index = False)