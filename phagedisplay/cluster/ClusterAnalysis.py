from Cluster import *

class ClusterAnalysis(Cluster):
	"""
	child class of Cluster.
	"""

	def elbowAnalysis(data):
	    """
	    find optimal k from data using elbow analysis.
	    """
	    initial = [vq.kmeans(data,i) for i in range(1,10)]
	    plt.plot([var for (cent,var) in initial])
	    plt.show()
	    
	def simplek(data):
	    """
	    super basic k approximation.
	    """
	    
	    k = np.sqrt(len(data)/2)
	    
	    return k

	def silhouette(data, n_clust):
	    """
	    return silhouette score
	    """
	    
	    calc = KMeans(n_clusters = n_clust).fit(data)
	    sil_score = silhouette_score(calc, calc.labels_, sample_size = len(data))
	    
	    return sil_score