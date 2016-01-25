import pyximport
pyximport.install()
import DistMatrix

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

class ClusterAnalysis():
    """
    analysis operations to determine best clustering.
    """

    def simplek(data):
        """
        super basic k approximation.
        """
        
        k = sqrt(len(cluster)/2)
        
        return k

    def silhouette(data, start, stop):
        """
        return silhouette score.
        Best score is closest to 1, worst score closest to -1, and scores near 0 indicates overlapping clusters.
        """
        
        s = []

        for n_clust in range(start, stop):
            calc = AgglomerativeClustering(n_clusters = n_clust, affinity = 'precomputed', linkage = 'average').fit(data.as_matrix()) 
         
                  
            labels = calc.labels_
            
            sil_score = silhouette_score(data.as_matrix(), labels, metric = 'precomputed')
            s.append(sil_score)
        
        plt.plot(s)
        plt.ylabel("Silhouette Score")
        plt.xlabel("k")
        plt.xlim((start, stop))
        plt.title("Cluster analysis for Agglomerative)

class Membership():
    """
    accessing cluster membership of given cluster in a pandas dataframe format.
    """
    
    def __init__(self, cluster):
        self._cluster = cluster.groupby('Cluster')
        
    def getAll(self):
        """
        return dictionary of all clusters and members.
        """
        
        return self._cluster.groups
    
    def getCluster(self, clust_num):
        """
        return specified cluster and members.
        """
        
        return self._cluster.get_group(clust_num)
    
    def getCount(self):
        """
        return count in each cluster.
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