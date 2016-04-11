"""
Main cluster processing class.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-04-11"

from phagedisplay import BaseProcessor
from phagedisplay import processors.cluster as c
    
class ClusterProcessor(BaseProcessor):

    def process(self,regession_out_dict,
                dist_function="weighted",
                dist_matrix="weight_matrices/blosum62.txt",
                cluster_method="DBSCAN"):

        self._dist_functions = {"hamming":c.dist_matrix.HammingDistMatrix,
                                "damerau":c.dist_matrix.DamerauDistMatrix,
                                "weighted":c.dist_matrix.WeightedDistMatrix}

        self._cluster_methods = {"DBSCAN":c.cluster.ClusterDB,
                                 "heirarchical":c.cluster.ClusterHCL}

        self._dist_matrix = self._dist_functions[dist_function]
        self._clusters = self._cluster_methods[cluster_method] 
     
        self._dist_matrix.create_data_vector(regression_out_dict)
        self._dist_matrix.calc_dist_matrix()
        self._dist_matrix.save_matrix()
 
        self._clusters = cluster.ClusterDB()
        self._clusters.generate_clusters(self._dist_matrix)
        self._clusters.write_output()

    @property
    def data(self):

        return self._C.cluster_labels
