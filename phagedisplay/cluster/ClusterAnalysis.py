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
	    X = squareform(data)
        linkage = hcl.average(X)

	    for n_clust in range(start, stop):
	    	calc = hcl.fcluster(linkage, self._factor, criterion = 'distance')
	     
	              
	        labels = calc.labels_
	        
	        sil_score = silhouette_score(data.as_matrix(), labels, metric = 'precomputed')
	        s.append(sil_score)
	    
	    plt.plot(s)
	    plt.ylabel("Silhouette Score")
	    plt.xlabel("t")
	    plt.xlim((start, stop))
	    plt.title("Cluster analysis for Agglomerative)