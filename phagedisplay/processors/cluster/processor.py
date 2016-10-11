"""
Main cluster processing class.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-04-11"

import os 
from .. import BaseProcessor
from . import dist_matrix, cluster 

class ClusterProcessor(BaseProcessor):

    def process(self,regression_out_dict,
                dist_function="weighted",
                weight_matrix_name="blosum62",
                cluster_method="DBSCAN"):

        self._dist_functions = {"hamming":dist_matrix.HammingDistMatrix,
                                "damerau":dist_matrix.DamerauDistMatrix,
                                "weighted":dist_matrix.WeightedDistMatrix}

        self._cluster_methods = {"DBSCAN":cluster.ClusterDB,
                                 "heirarchical":cluster.ClusterHCL}

        self._dist_matrix = self._dist_functions[dist_function](weight_matrix_name)
        self._clusters = self._cluster_methods[cluster_method] 
     
        self._dist_matrix.create_data_vector(regression_out_dict)
        self._dist_matrix.calc_dist_matrix()

        base_dir = self.getProperty("expt_name")
        self._clusters = cluster.ClusterDB(os.path.join(base_dir,"clusters"))
        self._clusters.generate_clusters(self._dist_matrix)
        self._clusters.write_output()

        self.saveFile(overwrite=True)


    @property
    def data(self):

        return self._clusters.cluster_labels
