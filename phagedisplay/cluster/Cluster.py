from sklearn.cluster import *
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform

import numpy as np


class Cluster():
    """
    Performs a given clustering algorithm on a distance matrix.
    args
        dist_matrix: distance matrix to perform clustering on
        cluster_alg: clustering algorithm to use on matrix
        factor: for setting a parameter on calculating clusters/number of predicted clusters.
        
    """
    
    def __init__(self, dist_matrix, cluster_alg, factor = None):
        self._dist_matrix = dist_matrix
        self._cluster_alg = cluster_alg
        self._factor = factor

    def cluster(self):
        """
        cluster using DBSCAN or hierarchal clustering. DBSCAN works best on larger datasets
        while hierarchal works best on smaller datasets.
        """

        X = self._dist_matrix

        # compute scipy hierarchal clustering
        if self._cluster_alg == 'hierarchal':
            condensed = squareform(X)
            linkage = hcl.average(condensed)
            clusters = hcl.fcluster(linkage, self._factor, criterion = 'distance')

        # compute DBSCAN clustering
        elif self._cluster_alg == 'DBSCAN':
            db = DBSCAN(eps = self._factor, metric = 'precomputed').fit(X.as_matrix() )
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