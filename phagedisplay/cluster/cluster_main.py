from DistMatrix import *
#from Cluster import *


#phagefile = sim_clusters #input("enter the file name: ")
#matrix = blosum62#input("enter the scoring matrix name: ")
#clustering = input("input clustering algorithm: ")

test = DistMatrix('hA5-0.0.txt', 'hamming', False).calcDistMatrix()
print(test)

#cluster tests
#Cluster(test, 'k_means', 3).silhouetteAnalysis()

