import numpy as np

def readFile(fileName):
	"""
	read matrix text file into dictionary for scoring.
	"""

	data = open(fileName).readlines()
	seq = data[0].strip('/n/r').split()

	matrix = {}

	for line in data[1:]:
		line = line.strip('/n/r').split()
		for j in range(1, len(line)):
			b = seq[j-1]
			matrix[(line[0], b)] = int(line[j])

	return matrix