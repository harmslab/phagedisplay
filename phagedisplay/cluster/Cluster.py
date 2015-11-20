from sklearn.cluster import *
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
import scipy.cluster.hierarchy as hcl

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


class Cluster():
    """
    Performs a given clustering algorithm on a distance matrix.

    args
        dist_matrix: distance matrix to perform clustering on
        cluster_alg: clustering algorithm to use on matrix
        factor: for setting a parameter on calculating clusters.
        num (optional): number of clusters predicted.
        
    """
    
    def __init__(self, dist_matrix, cluster_alg, factor, num = None):
        self._dist_matrix = dist_matrix 
        self._cluster_alg = cluster_alg
        self._factor = factor
        self._num = num
        
    def plotCluster(self, x, y):
        """
        plot cluster from algorithm.

        args
            x 
        """
        colors = plt.cm.Spectral(np.linspace(0, 1, len(self.cluster().clusters)))
        for k, col in zip(size, colors):
            
            #DBSCAN, use black for outlier
            if k == -1:
                col = 'k'
                
            cluster_centers = 
                
            members = labels == k
            plt.plot(X[members, 0], X[members, 1], 'o', markerfacecolor = col,
                     markeredgecolor = 'k', markersize = 14)
            
            plt.plot()
            plt.ylim(y)
            plt.xlim(x)

        plt.title('Estimated number of clusters: %d' % n_clusters_)
        plt.show()

    def cluster(self):
        """
        cluster and shit.
        """

        X = matrix

        # compute hierarchal clustering
        if self._cluster_alg == 'hierarchal':
            link = hcl.linkage(squareform(X)) # convert to condensed matrix before computing
            clusters = hcl.fcluster(link, self._factor, 'distance')

        # compute DBSCAN clustering
        elif self._cluster_alg == 'DBSCAN':
            db = DBSCAN(eps = self._factor, metric = 'precomputed').fit(X)
            cluster_centers = np.zeros_like(db.labels_, dtype=bool)
            cluster_centers[db.core_sample_indices_] = True
            clusters = db.labels_

        # compute mean shift clustering
        elif self._cluster_alg == 'mean_shift':
            bandwidth = estimate_bandwidth(X, quantile=self._factor)
            ms = MeanShift(bandwidth=bandwidth).fit(X)
            clusters = ms.labels_
            cluster_centers = ms.cluster_centers_

        # compute k-means clustering
        elif self.cluster_algo == 'k_means':
            kmeans = KMeans(n_clusters = self._num)
            clusters = kmeans.fit(X)
            cluster_centers = kmeans.cluster_centers_

        else:
            raise ValueError('algorithm unavailable')

        # get number of clusters
        n_clusters = len(np.unique(clusters)) - (1 if -1 in labels else 0)