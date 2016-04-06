__description__ = \
"""
"""
__author__ = "Hiranmayi Duvvuri"
__date__ = "2016-04-06"

import sklearn
from sklearn import cluster
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import pickle
import makeweblogos as logo

class Cluster:
    """
    Base class for clustering distance matrices (DistMatrix instances). 
    """

    def __init__(self):
        """
        Dummy init function.
        """

        pass

    def generate_clusters(self):
        """
        Dummy generate clusters function. This should be defined for each
        daughter class.
        """

        self.clust_labels = None

    def write_output(self,out_path):
        """
        Write out cluster output.
        """

        # save cluster size summary
        count = self.clust_labels.count()
        count.to_pickle('{}/summary.pickle'.format(out_path))
    
        # save each cluster in csv files
        for i in range(self.num_clusters):
            csv_file="{}/{}.csv".format(out_path,i)
            self.cluster_labels.get_group(i)["sequences"].to_csv(csv_file, index = False)
        
        # save weblogo of each cluster
        for i in range(self.num_clusters):
            logo.create_weblogo(self.cluster_labels.get_group(i)['sequences'].tolist(),
                                '{}/clust{}.pdf'.format(out_path, i))


class ClusterDB(Cluster):
    """
    Generate clsuters by dbscan.
    """

    def __init__(self,min_samples,epsilon=None,metric="precomputed",leaf_size=300,
                 algorithm="ball_tree",epsilon_size_cutoff=0.95):
        """
        Initialize clustering options.  The options nature of each option can
        be found in the scikitlearn dbscan documention.
        """

        self.min_samples = min_samples
        self.epsilon = epsilon
        self.metric = metric
        self.leaf_size = leaf_size
        self.algorithm = algorithm
        self.epsilon_size_cutoff = epsilon_size_cutoff
   
    def generate_clusters(self,D):
        """
        Generate clusters by dbscan.  Takes a DistMatrix instance D as input.
        """

        # If the user does not specify epsilon, estimate it.  
        if self.epsilon == None:
            self.epsilon = self._estimate_epsilon(D.dist_frame)

        # Initialize a cluster.DBSCAN instance
        db = cluster.DBSCAN(eps=self.epsilon,
                            min_samples=self.min_samples,
                            metric=self.metric, 
                            leaf_size=self.leaf_size,
                            algorithm=self.algorithm)

        # Generate clusters
        db.fit(D.dist_frame.as_matrix())

        # Grab core samples
        self.core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        self.core_samples_mask[db.core_sample_indices_] = True
        self.clusters = db.labels_

        # Basic cluster stats
        self.num_clusters = n_clusters = len(np.unique(self.clusters)) - (1 if -1 in clusters else 0)
        self.cluster_labels = pd.DataFrame({'sequences' : D.dist_frame.index, 
                                            'cluster' : self.clusters})

    def _estimate_epsilon(self,D)
        """
        Choose the value of epsilon that maximizes both number of clusters 
        (N/Nmax > epsilon_size_cutoff) and maximizes the size of the noise
        cluster.
        """

        epsilon_list = []
        num_clust_list = []
        noise_list = []

        # Go through a large number of values of epsilon 
        for i in range(2,int(max(D.dist_frame.max()))):

            # generate clusters at this value of epsilon
            self.epsilon = i
            self.generate_clusters(D.dist_frame)

            # record the epsilon, number of clusters, and size of the noise cluster
            epsilon_list.append(i)
            num_clust_list.append(self.num_clusters)
            noise_list.append(len(self.cluster_labels[(self.cluster_labels['cluster'] == -1)].index))

            # if we created more than just one giant cluster, write out the membership
            if self.num_clusters > 1:
                count = self.cluster_labels.count()
                count.to_pickle("junk_{:d}.pickle".format(i))

        # Normalize the number of clusters to the largest number seen
        clust_thresh = np.array(num_clust_list)/max(num_clust_list)

        # Get indices of each epsilon where the number of clusters is above
        # epsilon_size_cutoff.
        indices = np.where(clust_thresh[clust_thresh > self.epsilon_size_cutoff])

        # Now find values of epsilon that maximize the size of the noise cluster
        max_noise = max(noise_list)
        eps = [epsilon_list[i] for i in indices if noise_list[i] == max_noise]

        return eps[0]

class ClusterHCL(Cluster):
    """
    Cluster by heirarchical clustering. 
    """

    def __init__(self,factor,criterion="maxclust"):
        """
        Perform heirarchical clustering.  See the scipy heirarchical clustering
        documentation.
        """
        
        self.factor = factor
        self.criterion = criterion
    
    def generate_clusters(self,D):
        """
        Generate clusters.  Takes a DistMatrix instance D. 
        """

        condensed = squareform(D.dist_frame)
        linkage = hcl.average(condensed)
        self.clusters = hcl.fcluster(linkage,self.factor,criterion=self.criterion)

        self.num_clusters = n_clusters = len(np.unique(self.clusters)) - (1 if -1 in clusters else 0)
        self.cluster_labels = pd.DataFrame({'sequences' : D.dist_frame.index, 
                                            'cluster' : self.clusters})
    

class EpsAnalysis():
    """
    Analyze epsilon.
    """   
 
    def graphs(self):
        """
        return epsilon and noise vs clusters graphs.
        """

        with PdfPages('Cluster_test/all_epsilons/eps_noise.pdf') as pdf:

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
    
        
class ClusterProcessor(BaseProcessor):
    """
    """

    def process(self):
        """
        """

        pass

    @property
    def data(self):
        """
        """

        return None
