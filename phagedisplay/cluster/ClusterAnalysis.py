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