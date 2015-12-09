import numpy as np

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