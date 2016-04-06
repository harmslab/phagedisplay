__description__ = \
"""
"""
__author__ = "Hiranmayi Duvvuri"
__date__ = "2016-04-06"

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
    make histogram plot of frequency of sequence distance scores in the
    given distance matrix.
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
    
    def __init__(self, dist_matrix, cluster_alg, factor=None,num=None):

        self._dist_matrix = dist_matrix
        self._cluster_alg = cluster_alg
        self._factor = factor
        self._num = num

    def cluster(self):
        """
        Cluster a distance matrix.  

        DBSCAN - 
            min_samples: the minimum number of samples for the point to be
                         considered to be a central cluster point. A higher
                         number will give a smaller number of large clusters
                         with a large noise (-1) cluster while a lower number
                         will give a large number of small clusters with a
                         very small or no noise cluster.
            epsilon:     use of freqHist can help determine a factor to use,
                         max distance radius around a point to grab cluster
                         points from. 
        """

        X = self._dist_matrix
        
        # compute scipy hierarchal clustering
        if self._cluster_alg == 'hierarchal':
            condensed = squareform(X)
            linkage = hcl.average(condensed)
            clusters = hcl.fcluster(linkage, self._factor, criterion = 'maxclust')

        # compute DBSCAN clustering
        elif self._cluster_alg == 'DBSCAN':

            db = DBSCAN(eps = self._factor,
                        min_samples = self._num,
                        metric = 'precomputed', 
                        leaf_size = 300,
                        algorithm = 'ball_tree').fit(X.as_matrix())

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
    """
    XX What does this do?
    """    
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
    
    def epsilon(self, dist_matrix, file, cutoff = 0.95):
        """
        return epsilon choice based on threshold N/Nmax > cutoff and the maximum noise.
        """
    
        self.graphs()
    
        for i in range(2, int(max(dist_matrix.max()))):

            clustering = Cluster(dist_matrix,
                                 'DBSCAN',
                                 factor = i,
                                 num = 7).cluster()

            length = self.num_clust(clustering)
            self._epsilon.append(i)
            self._clusters.append(length)
        
            outliers = len(clustering[(clustering['Cluster'] == -1)].index)
            self._noise.append(outliers)

            self._all_combo.append((i, length, outliers))
        
            if length > 1:
                # cluster summary for each epsilon
                count = Membership(clustering).get_count()
                count.to_pickle('Cluster_test/all_epsilons/eps{}.pickle'.format(i))

        data = np.array(self._all_combo)

        # Normalization 
        max_clust = max(data[:,1])

        # Normalize the number of clusters
        clust_thresh = data[:,1]/max_clust

        # Get indices in array that satisfy condition 1
        indices = np.where(clust_thresh[clust_thresh > cutoff])
        max_noise = max(data[indices, 2])

        # Get indices that satisfy condition 2
        eps = [self._all_combo[i][0]
               for i in indices if self._all_combo[i][2] == max_noise]

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
    count.to_pickle('{}/summary.pickle'.format(summary_file))
    
    # save each cluster in csv files
    for i in range(clust_num):
        clust_mem.to_csv(i, '{}/{}.csv'.format(clust_loc, i))
     
    # save weblogo of each cluster
    for i in range(clust_num):
        logo.create_weblogo(clust_mem.get_cluster(i)['Sequences'].tolist(), '{}/clust{}.pdf'.format(logo_loc, i))


class ClusterProcessor(BaseProcessor):
    """
    """

    def process(self):
        """
        Do the regression. 
            count_dict is a set of sequences with counts over rounds as values
            cg_width is the number of standard deviations at which counts in a
                     given round are considered different
            minium_times_seen is the number of times a particular sequence must
                              be seen across all rounds if it is to be included
                              in the fit.
            human_out_file human-readable output file
            global_regression (True/False) whether or not to run global     
                              regression rather than simply fitting each pattern
                              individually.
        """

        self.count_dict = count_dict
        self.cg_width = cg_width
        self.minimum_times_seen = minimum_times_seen

        # Create patterns
        self._find_patterns()

        # Set up fit
        fit_pattern = np.zeros((len(self.patterns),self.patterns[0].pattern_length),dtype=float)
        degeneracy = np.zeros((len(self.patterns)),dtype=int)
        for i, c in enumerate(self.patterns):
            if len(c.degen_sequence_set) != 0:
                mean_pattern, degen = c.calcMeanPattern()
                fit_pattern[i,:] = mean_pattern
                degeneracy[i] = degen

        # Do actual fit
        m = FitModel(fit_pattern,degeneracy,self.rounds,self._log_file)
        
        if global_regression:
            m.runRegression()

            # Grab values from fit
            theta, K, calc_values = m.returnParam()

        else:
            theta = np.array([0.0 for i in range(len(self.patterns))])
            K     = np.array([0.0 for i in range(len(self.patterns))])
            calc_values = np.array([[0.0 for j in range(max(self.rounds))]
                                    for i in range(len(self.patterns))])

        # create output
        self._out_dict = {}
        for i, c in enumerate(self.patterns):
            if len(c.degen_sequence_set) != 0:
                for j in range(len(c.degen_sequence_set)):
                    for s in c.degen_sequence_set[j].sequences:
                        self._out_dict[s] = (np.exp(K[i]),
                                             np.exp(theta[i]/degeneracy[i]),
                                             np.exp(m.param_guess[i+m.num_patterns]),
                                             np.exp(m.param_guess[i]/degeneracy[i]))

        new_out = [] 
        for k in self._out_dict.keys():
            new_out.append((tuple(self._out_dict[k]),k))
        new_out.sort(reverse=True)

        base_dir = self.getProperty("expt_name")
        f = open(os.path.join(base_dir,human_out_file),"w")
        f.write("{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}\n".format(" ",
                                                                      "sequence",
                                                                      "K_global",
                                                                      "theta_global",
                                                                      "K_indep",
                                                                      "theta_indep"))

        for i, o in enumerate(new_out):
            f.write("{:20d}{:>20s}{:20.10e}{:20.10e}{:20.10e}{:20.10e}\n".format(i,
                                                                                 o[1],
                                                                                 o[0][0],
                                                                                 o[0][1],
                                                                                 o[0][2],
                                                                                 o[0][3]))
        f.close()


    def _find_patterns(self):
        """
        Load the data into a set of degenerate patterns.
        """
        
        all_sequences = list(self.count_dict.keys())

        # Create a list, (degen_sequence_list) that has a list of all sequences that
        # share exactly the same count pattern.  Each unique pattern/sequence set is
        # stored in an instance of the DegenerateSequences class.
        self.degen_sequence_dict = {}
        self.rounds = [i for i, v in enumerate(self.count_dict[all_sequences[0]])
                       if v != None]
        num_rounds = len(self.rounds)
        for s in all_sequences:

            # Get rid of missing data
            pattern = tuple([v for v in self.count_dict[s] if v != None])

            if len(pattern) != num_rounds:
                err = "All sequences must have same rounds observed."
                raise ValueError(err)
           
            # If the sequence was *only* seen in the starting library, ignore it 
            if sum(pattern[1:]) == 0:
                continue
    
            if sum(pattern) < self.minimum_times_seen:
                continue
            try:
                self.degen_sequence_dict[pattern].append(s)
            except KeyError:
                self.degen_sequence_dict[pattern] = DegenerateSequences(pattern,s)

        degen_sequence_list = []
        for k in self.degen_sequence_dict.keys():
            degen_sequence_list.append(copy.copy(self.degen_sequence_dict[k]))

        num_unique_patterns = len(degen_sequence_list)
        if num_unique_patterns == 0:
            err = "No sequences were seen enough times to be included."
            raise ValueError(err)

        self._logger("Found {:d} unique patterns for {:d} clones".format(num_unique_patterns,
                                                                         len(all_sequences)))

        # Create a numpy array of all unique patterns
        unique_patterns = np.zeros((num_unique_patterns,num_rounds),dtype=int)
        for i in range(num_unique_patterns):
            unique_patterns[i,:] = np.array(degen_sequence_list[i].pattern)

        # Coarse grain these unique patterns into similar patterns
        cg = CoarseGrainer(np.max(unique_patterns),self.cg_width)
        similar_pattern_dict = {}
        for i in range(num_unique_patterns):

            key = tuple(cg.coarseGrain(unique_patterns[i,:]))
            try:
                similar_pattern_dict[key].append(degen_sequence_list[i])
            except KeyError:
                similar_pattern_dict[key] = [degen_sequence_list[i]]

        # Create a final set of clusters (made of SequenceCluster instance)
        self.patterns = [DegenerateSequenceCluster()
                         for i in range(len(similar_pattern_dict))]
        for i, key in enumerate(similar_pattern_dict.keys()):
            for cluster in similar_pattern_dict[key]:
                self.patterns[i].addSequenceSet(cluster)

        self._logger("Collapsed to {:d} similar patterns.".format(len(self.patterns)))

    @property
    def data(self):
        """
        Return a dictionary of sequences keyed to fit affiniites and initial counts.
        """

        return self._out_dict
