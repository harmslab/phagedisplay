"""
Data structures for clustering and visualizing sequences.
"""
__author__ = "Hiranmayi Duvvuri, Michael J. Harms"
__date__ = "2016-04-06"

import sklearn
from sklearn import cluster
import scipy.cluster.hierarchy as hcl
from scipy.spatial.distance import squareform

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import pickle, sys, os
import weblogolib, corebio

def create_weblogo(seq_list,out_file,alphabet="amino"):
    """ 
    Create a pdf sequence logo given a list of sequences.
    """

    # Alphabet and color scheme
    if alphabet == "amino":
        wl_alphabet = corebio.seq.unambiguous_protein_alphabet
        wl_color_scheme = weblogolib.colorscheme.hydrophobicity
    elif alphabet == "dna":
        wl_alphabet = corebio.seq.unambiguous_dna_alphabet          
        wl_color_scheme = weblogolib.colorscheme.nucleotide
 
    # load data and create a logo
    logo_sequences = corebio.seq.SeqList(seq_list,alphabet=wl_alphabet)
    logo_data = weblogolib.LogoData.from_seqs(logo_sequences)

    # create format
    logo_format = weblogolib.LogoFormat(logo_data)
    logo_format.color_scheme = wl_color_scheme

    # dump out as a pdf.    
    logo = weblogolib.pdf_formatter(logo_data, logo_format)
    f = open(out_file,'w')
    f.buffer.write(logo)
    f.close()

class Cluster:
    """
    Base class for clustering distance matrices (DistMatrix instances). 
    """

    def __init__(self,out_path,verbose=False):
        """
        Dummy init function.
        """

        self.out_path = out_path
        if not os.path.exists(out_path):
            os.mkdir(out_path) 

        self.verbose = verbose

    def generate_clusters(self):
        """
        Dummy generate clusters function. This should be defined for each
        daughter class.
        """

        self.cluster_labels = None

    def write_output(self,alphabet="amino"):
        """
        Write out cluster output.
        """

        print("Creating output."); sys.stdout.flush()

        g = self.cluster_labels.groupby("cluster")
        zero_pad_size = len("{:d}".format(self.num_clusters))

        # save cluster size summary
        count = g.count() 
        count.to_csv(os.path.join(self.out_path,"summary_count.txt"))

        for i in range(self.num_clusters):

            # save each cluster in csv files
            num = "{:d}".format(i).zfill(zero_pad_size)

            csv_string = "{}_cluster.csv".format(num)
            csv_file=os.path.join(self.out_path,csv_string)
            g.get_group(i)["sequences"].to_csv(csv_file,index=False)
        
            # save weblogo of each cluster
            pdf_string = "{}_cluster.pdf".format(num)
            pdf_file = os.path.join(self.out_path,pdf_string)

            create_weblogo(g.get_group(i)['sequences'].tolist(),
                           pdf_file,
                           alphabet)

class ClusterDB(Cluster):
    """
    Generate clsuters by dbscan.
    """

    def __init__(self,out_path,min_samples=5,epsilon=None,metric="precomputed",leaf_size=300,
                 algorithm="ball_tree",epsilon_size_cutoff=0.95,verbose=False):
        """
        Initialize clustering options.  The options nature of each option can
        be found in the scikitlearn dbscan documention.
        """
        
        super(self.__class__,self).__init__(out_path,verbose)

        self.out_path = out_path
        self.min_samples = min_samples
        self.epsilon = epsilon
        self.metric = metric
        self.leaf_size = leaf_size
        self.algorithm = algorithm
        self.epsilon_size_cutoff = epsilon_size_cutoff
        self.verbose = verbose
   
    def generate_clusters(self,D):
        """
        Generate clusters by dbscan.  Takes a DistMatrix instance D as input.
        """
        # If the user does not specify epsilon, estimate it.  
        if self.epsilon == None:
            self.epsilon = self._estimate_epsilon(D)

        # Initialize a cluster.DBSCAN instance
        db = cluster.DBSCAN(eps=self.epsilon,
                            min_samples=self.min_samples,
                            metric=self.metric, 
                            leaf_size=self.leaf_size,
                            algorithm=self.algorithm)
      
        db.fit(D.dist_frame.as_matrix())

        # Grab core samples
        self.core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        self.core_samples_mask[db.core_sample_indices_] = True
        self.clusters = db.labels_

        # Basic cluster stats
        self.num_clusters = len(np.unique(self.clusters)) - (1 if -1 in self.clusters else 0)
        self.cluster_labels = pd.DataFrame({'sequences' : D.dist_frame.index, 
                                            'cluster' : self.clusters})

    def write_output(self,alphabet="amino"):
        """
        Write out cluster output.
        """

        super(self.__class__,self).write_output(alphabet)

        f = open(os.path.join(self.out_path,"cluster_stats.txt"),"w")
        f.write("Clustered by dbscan\n")
        f.write("min_samples: {}\n".format(self.min_samples))
        f.write("epsilon: {}\n".format(self.epsilon))
        f.write("metric: {}\n".format(self.metric))
        f.write("leaf_size: {}\n".format(self.leaf_size))
        f.write("algorithm: {}\n".format(self.algorithm))
        f.write("epsilon_size_cutoff: {}\n".format(self.epsilon_size_cutoff))
        f.write("num_clusters: {}\n".format(self.num_clusters))
        f.close()

    def _estimate_epsilon(self,D):
        """
        Choose the value of epsilon that maximizes both number of clusters 
        (N/Nmax > epsilon_size_cutoff) and maximizes the size of the noise
        cluster.
        """
            
        print("Optimizing epsilon."); sys.stdout.flush()

        epsilon_list = []
        num_clust_list = []
        noise_list = []

        # Go through a large number of values of epsilon 
        for i in np.arange(0,np.max(D.dist_matrix),0.1):

            # generate clusters at this value of epsilon
            self.epsilon = i

            # This check is because dbscan throws an error if epsilon is too small...
            try:
                self.generate_clusters(D)
            except ValueError:
                continue

            # record the epsilon, number of clusters, and size of the noise cluster
            epsilon_list.append(i)
            num_clust_list.append(self.num_clusters)
            noise_list.append(len(self.cluster_labels[(self.cluster_labels['cluster'] == -1)].index))

            # spit out epsilon optimization if being verbose
            if self.verbose:
                print(epsilon_list[-1],num_clust_list[-1],noise_list[-1])
                sys.stdout.flush()
            
                if self.num_clusters > 1:
                    count = self.cluster_labels.groupby("cluster").count()
                    count.to_pickle(os.path.join(self.out_path,"episilon_{:.2e}.pickle".format(i)))

        # If no clusters were found for *any* epsilon, complain
        if len(num_clust_list) < 1:
            err = "No clusters found for any epsilon.  Data set has too few sequences?\n"
            raise ValueError(err)

        # Normalize the number of clusters to the largest number seen
        clust_thresh = np.array(num_clust_list)/max(num_clust_list)

        # Get indices of each epsilon where the number of clusters is above
        # epsilon_size_cutoff.
        indices = np.where(clust_thresh > self.epsilon_size_cutoff)

        # Now find values of epsilon that maximize the size of the noise cluster
        max_noise = max([noise_list[i] for i in indices[0]])
        eps = [epsilon_list[i] for i in indices[0] if noise_list[i] == max_noise]
  
        # return the smallest epsilon compatible with this.
        return eps[0]

class ClusterHCL(Cluster):
    """
    Cluster by heirarchical clustering. 
    """

    def __init__(self,out_path,factor,criterion="maxclust",verbose=False):
        """
        Perform heirarchical clustering.  See the scipy heirarchical clustering
        documentation.
        """

        super(self.__class__,self).__init__(out_path,verbose)
        
        self.out_path = out_path
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
                                            'cluster'   : self.clusters})
    
