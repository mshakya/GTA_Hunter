"""
File name: Weight.py
Date created: 01/12/2016
Date last modified: 04/07/2016
Python version: 3.5.1
Description: Weighting system used to 
	adjust SVM classifcation of genes
	during training.
"""

###############
### IMPORTS ###
###############
from Profile import Profile
from Loader import Loader
import numpy as np
import matplotlib.pyplot as plt
import argparse

############
### CODE ###
############
def flatten(x):
		"""List flattening function taken
			from stackoverflow.com 
		"""
		result = []
		for el in x:
			if hasattr(el, "__iter__") and not isinstance(el, str):
				result.extend(flatten(el))
			else:
				result.append(el)
		return result

class Weight:
	"""Used to create weights for training set"""

	def __init__(self, profiles, pairwiseDict):
		""" Init for Weight. Takes Loader of the
			profiles.
			Input:
				profiles: Loader of profiles
				pairwiseDict: a 2D dictionary of
					gene names
			Returns:
				None.
		"""
		self.profiles = profiles
		self.pairwiseDict = pairwiseDict

	@staticmethod
	def load(filename):
		"""Loads .dist file from RAxML output
			and creates a dictionary of 
			pairwise distances, which it 
			returns.
		"""
		pairwiseDict = {}
		for line in open(filename, 'r'):
			parts = line.split()
			if parts[0] in pairwiseDict:
				pairwiseDict[parts[0]][parts[1]] = float(parts[2])
			else:
				pairwiseDict[parts[0]] = {parts[1]: float(parts[2])}
		return pairwiseDict

	def cluster(self, clusterType, cutoff):
		"""Clusters profiles after distance has
			been calculated. Hierarchical clustering
			can be done using 'farthest' or 'nearest'
			neighbors (type).
		"""
		if clusterType == "nearest":
			multiplier = 1
		elif clusterType == "farthest":
			multiplier = -1
		else:
			print("Cluster type",clusterType,"not recongnized.")
			return None
		# Build list of clusters
		clusters = [[profile] for profile in self.profiles]
		# Loop until all clusters are merged
		while len(clusters) > 1:
			# Init
			minClustDist = float("inf")
			minClust1 = None
			minClust2 = None
			# Compare cluster pairs
			for cluster1 in clusters:
				for cluster2 in clusters:
					# Compute distance if different
					if cluster1 != cluster2:
						# Loop over profiles in clusters
						bestDist = multiplier * float("inf")
						bestProf1 = None
						bestProf2 = None
						for profile1 in flatten(cluster1):
							for profile2 in flatten(cluster2):
								# Lookup profile distance
								if profile1.name in self.pairwiseDict and \
									profile2.name in self.pairwiseDict[profile1.name]:
									currentDist = self.pairwiseDict[profile1.name][profile2.name]
								else: # Other order
									currentDist = self.pairwiseDict[profile2.name][profile1.name]
								# Update best match (farthest or nearest profiles between clusters)
								if multiplier > 0 and currentDist < bestDist: # min/nearest
									bestDist = currentDist
									bestProf1 = profile1
									bestProf2 = profile2
								elif multiplier < 0 and currentDist > bestDist: # max/farthest
									bestDist = currentDist
									bestProf1 = profile1
									bestProf2 = profile2
						# Once distance for profiles found, get closest cluster
						if bestDist < minClustDist:
							minClustDist = bestDist
							minClust1 = cluster1
							minClust2 = cluster2
			#Check distance cutoff
			if minClustDist > cutoff:
				#Exceeded cutoff, end clustering
				return clusters
			else:
				# Merge clusters and remove others
				clusters = clusters + [[minClust1, minClust2]]
				clusters.remove(minClust1)
				clusters.remove(minClust2)

		return clusters

	def weight(self, clusters):
		"""Weights each profile based on
			how they are clustered together
		"""
		for cluster in clusters:
			profiles = flatten(cluster)
			cluster_size = len(profiles)
			for profile in profiles:
				profile.weight = 1.0 / cluster_size

	""" Visualization and Analysis Functions """
	def visualize_clusters(self, clusters):
		for cluster in clusters:
			self.visualize_helper(cluster, True, ' ')
			print("-*")
			
	def visualize_helper(self, cluster, isFirst, indent):
		angle = '/' if isFirst else '\\'
		if len(cluster) ==1:
			print(indent+' '+angle+'-'+str(cluster[0]))
		else:
			self.visualize_helper(cluster[0], True, indent + ('   ' if isFirst else ' | '))
			print(indent+' '+angle+'-*')
			self.visualize_helper(cluster[1], False, indent + (' | ' if isFirst else '   '))

	def plot_hist(self):
		vals = []
		for key1 in pairwiseDict:
			for key2 in pairwiseDict[key1]:
				vals.append(pairwiseDict[key1][key2])

		plt.hist(vals)
		plt.title("Pairwise Distance Histogram")
		plt.show()

	def print_quartiles(self):
		vals = []
		for key1 in pairwiseDict:
			for key2 in pairwiseDict[key1]:
				vals.append(pairwiseDict[key1][key2])

		vals.sort()

		l = len(vals)

		print("Q1: %f - %f" % (vals[0], vals[int(l/4)]))
		print("Q2: %f - %f" % (vals[int(l/4)], vals[int(l/2)]))
		print("Q3: %f - %f" % (vals[int(l/2)], vals[int(3*l/4)]))
		print("Q4: %f - %f" % (vals[int(3*l/4)], vals[l-1]))

	def print_under_threshold(self, cutoff):
		vals = []
		for key1 in pairwiseDict:
			for key2 in pairwiseDict[key1]:
				vals.append(pairwiseDict[key1][key2])
		vals2 = [i for i in vals if i < cutoff]
		print("%d/%d under threshold of %f" % (len(vals2), len(vals), cutoff))

	def count_clusters(self, clusters):
		clustered = 0
		cl = 0
		for cluster in clusters:
			if len(cluster) > 1:
				clustered += len(flatten(cluster))
				cl += 1
		print("%d genes have been clustered, out of %d total, over %d clusters." % (clustered, len(self.profiles), cl))

# Command-line driver
if __name__ == '__main__':
	# Define arg parser
	parser = argparse.ArgumentParser(description="Gene weighting.")

	parser.add_argument("-p", "--profile_path", type=str, nargs=1,
		dest="profile_path", required=True,
		help="The .faa or .fna file used in calculating pairwise distances")
	parser.add_argument("-w", "--pairwise_path", type=str, nargs=1,
		dest="pairwise_path", required=True,
		help="The .dist file of pairwise distances created by RAxML")
	parser.add_argument("-t", "--cluter_type", type=str, nargs=1,
		dest="cluster_type", required=False, default=["farthest"],
		help="Specify 'farthest' or 'nearest' neighbors clustering.")
	parser.add_argument("-d", "--cutoff_distance", type=float, nargs=1,
		dest="cutoff", required=False, default=[0.05],
		help="Specify the cutoff distance for clustering.")
	parser.add_argument("-v", "--visualize", action="store_true",
		dest="v", required=False,
		help="Shows visualization of the clustered profiles.")
	parser.add_argument("-i", "--histogram", action="store_true",
		dest="i", required=False,
		help="Shows histogram of pairwise distances.")
	parser.add_argument("-q", "--quartiles", action="store_true",
		dest="q", required=False,
		help="Shows quartiles of pairwise distances.")
	parser.add_argument("-r", "--threshold", type=float, nargs=1,
		dest="r", required=False,
		help="Show number of pairwise distances that fall below the given threshold.")
	parser.add_argument("-c", "--count", action="store_true",
		dest="c", required=False,
		help="Shows how many profiles clustered into how many clusters.")

	args = parser.parse_args()
	
	# Load in profiles and pairwise dictionary
	profiles = Loader.load(args.profile_path[0])
	pairwiseDict = Weight.load(args.pairwise_path[0])
	# Init Weight object
	weight = Weight(profiles, pairwiseDict)
	# Create clusters
	clusters = weight.cluster(args.cluster_type[0], args.cutoff[0])
	# Assign weights to profiles
	weight.weight(clusters)
	# Visualize
	if args.v:
		weight.visualize_clusters(clusters)
	# Histogram
	if args.i:
		weight.plot_hist()
	# Quartiles
	if args.q:
		weight.print_quartiles()
	# Threshold
	if args.r:
		weight.print_under_threshold(args.r[0])
	# Count
	if args.c:
		weight.count_clusters(clusters)


	
	