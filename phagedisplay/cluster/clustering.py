from sklearn.cluster import DBSCAN
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import pickle
import makeweblogos as logo

def histogram(dist_matrix):
        """
        make histogram plot of frequency of sequence distance scores in the given distance matrix.
        """
        
        plt.figure();
        dist_matrix.plot(kind='hist', legend=False, orientation='horizontal')

class Cluster():
    """
    Performs a given clustering algorithm on a distance matrix.
    args
        dist_matrix: distance matrix to perform clustering on
        cluster_alg: clustering algorithm to use on matrix
        factor: for setting a parameter on calculating clusters.
        num (optional): number of clusters predicted.
        
    """
    
    def __init__(self, dist_matrix, cluster_alg, factor = None, num = None):
        self._dist_matrix = dist_matrix
        self._cluster_alg = cluster_alg
        self._factor = factor
        self._num = num

    def cluster(self):
        """
        clustering to perform on distance matrix.
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
            clusters = hcl.fcluster(linkage, self._factor, criterion = 'maxclust')

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
    
    
class EpsAnalysis():
    
    def __init__(self):
        self._epsilon = []
        self._clusters = []
        self._noise = []
        self._all_combo = []
        
    def num_clust(self, matrix_clust):
        """
        get number of clusters.
        """
        clusters = matrix_clust['Cluster'].tolist()
        num = len(np.unique(clusters)) - (1 if -1 in clusters else 0)
        
        return num
            
    def graphs(self, eps_file):
        """
        return epsilon and noise vs clusters graphs.
        """
        with PdfPages(eps_file) as pdf:
            plt.subplot(2, 1, 1)
            plt.plot(self._epsilon, self._clusters)
            plt.ylabel("clusters")
            plt.xlabel("epsilon")
            pdf.savefig()
            plt.close()

            plt.subplot(2, 1, 2)
            plt.plot(self._noise, self._clusters)
            plt.ylabel("clusters")
            plt.xlabel("noise")
            pdf.savefig()
            plt.close()
    
    def epsilon(self, dist_matrix, eps_file, eps_summary):
        """
        return epsilon choice based on threshold N/Nmax > 0.95 and the maximum noise.
        """
    
        self.graphs(eps_file)
    
        for i in range(2, int(max(dist_matrix.max()))):
            clustering = Cluster(dist_matrix, 'DBSCAN', factor = i, num = 7).cluster()
            length = self.num_clust(clustering)
            self._epsilon.append(i)
            self._clusters.append(length)
        
            outliers = len(clustering[(clustering['Cluster'] == -1)].index)
            self._noise.append(outliers)

            self._all_combo.append((i, length, outliers))
        
            if length > 1:
                # cluster summary for each epsilon
                count = Membership(clustering).get_count()
                count.to_pickle('{}/eps{}.pkl'.format(eps_summary, i))

        data = np.array(self._all_combo)

        # Normalization 
        max_clust = max(data[:,1])

        # Normalize the number of clusters
        clust_thresh = data[:,1]/max_clust

        # Get indices in array that satisfy condition 1
        indices = np.where(clust_thresh[clust_thresh > 0.95])
        max_noise = max(data[indices, 2])

        # Get indices that satisfy condition 2
        eps = [self._all_combo[i][0] for i in indices if self._all_combo[i][2] == max_noise]

        return eps[0]
        

class Membership():
    """
    accessing cluster membership of given cluster in a pandas dataframe format.
    """
    
    def __init__(self, cluster):
        self._cluster = cluster.groupby('Cluster')
        
    def get_all(self):
        """
        return dictionary of all clusters and members.
        """
        
        return self._cluster.groups
    
    def get_cluster(self, clust_num):
        """
        return specified cluster and members.
        """
        
        return self._cluster.get_group(clust_num)
    
    def get_count(self):
        """
        return count in each cluster.
        """
        
        return self._cluster.count()
    
    def get_overlap(self, x, y, clust2):
        """
        compare two clusters to see amount of overlap.
        """
        
        overlap = []

        for i in range(x):
            for j in range(y):
                a = self.get_cluster(i)
                b = clust2.get_cluster(j)
                percent_overlap = round(((len(np.intersect1d(a.as_matrix(['Sequences']), b.as_matrix(['Sequences'])))
                                         /len(a))*100), 2)

                # returns unique values in both arrays.
                overlap.append((i, j, len(a), len(b), percent_overlap))

        col = ['cluster a', 'cluster b', 'length a', 'length b', '% overlap']
        
        return  pd.DataFrame(overlap, columns = col)
    
    def to_csv(self, cluster, file):
        """
        save a cluster to a csv file.
        """
        
        self.get_cluster(cluster)['Sequences'].to_csv(file, index = False)
        
        
def cluster_auto(matrix, eps_file, eps_summary, summary_file, clust_loc, logo_loc):
    """
    clusters from distance matrix, cluster data output.
    """
    
    clust_test = Cluster(matrix, 'DBSCAN', factor = EpsAnalysis().epsilon(matrix, eps_file, eps_summary), num = 7).cluster()
    clust_mem = Membership(clust_test)
    
    clust_num = EpsAnalysis().num_clust(clust_test)
    
    # save cluster size summary
    count = clust_mem.get_count()
    count.to_pickle('{}/summary.pkl'.format(summary_file))
    
    # save each cluster in csv files
    for i in range(clust_num):
        clust_mem.to_csv(i, '{}/{}.csv'.format(clust_loc, i))
     
    # save weblogo of each cluster
    for i in range(clust_num):
        logo.create_weblogo(clust_mem.get_cluster(i)['Sequences'].tolist(), '{}/clust{}.pdf'.format(logo_loc, i))