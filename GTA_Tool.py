"""
File name: GTA_Tool.py
Date created: 10/22/2015
Date last modified: 04/07/2016
Python version: 3.5.1
Description: Combines Loader.py,
	Profile.py, SVM.py, Features.py,
	Weight.py, and Filter.py to be able
	to characterize and classify two different
	types of genes (GTA and phage).
"""

###############
### IMPORTS ###
###############
from Loader import Loader
from Weight import Weight
from Feature import Feature
from SVM import SVM
import argparse
import numpy as np
import sys

############
### META ###
############
__author__ = "Taylor Neely"
__version__ = "1.0.0"
__email__ = "tneely@dartmouth.edu"

#################
### CONSTANTS ###
#################
NREPS = 10 # for xval in SVM
PSE_WEIGHT = 0.05 # for pseaac feature

############
### CODE ###
############
def get_args():
	# Init parser
	parser = argparse.ArgumentParser(description="Gene Classification Using SVM.")

	### Define Args ###

	# Main
	parser.add_argument("-g", "--GTA", type=str, nargs=1,
		dest="gta", required=True,
		help="The .faa or .fna training file for GTA genes.")
	parser.add_argument("-v", "--virus", type=str, nargs=1,
		dest="virus", required=True,
		help="The .faa or .fna training file for viral genes.")
	parser.add_argument("-q", "--queries", type=str, nargs=1,
		dest="queries", required=False,
		help="The .faa or .fna query file to be classified.")
	parser.add_argument("-k", "--kmer", type=int, nargs=1,
		dest="kmer", required=False, default=[4],
		help="The kmer size needed for feature generation (default=4).")
	parser.add_argument("-p", "--pseaac", nargs="?",
		dest="pseaac", required=False, const=3, default=None,
		help="Expand feature set to include pseudo amino acid composition. Specify lamba (default=3) and weight (default = 0.05).")
	parser.add_argument("-m", "--min", action="store_true",
		dest="mini", required=False,
		help="Print bare minimum results.")

	# Weight
	parser.add_argument("-w", "--weight", type=str, nargs=2,
		dest="weight", required=False,
		help="Allows for weighting of training set. Will need to specify the two pairwise distance\
			files needed for weighting (GTA first, then virus).")
	parser.add_argument("-t", "--cluster_type", type=str, nargs=1,
		dest="cluster_type", required=False, default=['farthest'],
		help="Specify 'farthest' or 'nearest' neighbors clustering (default='farthest').")
	parser.add_argument("-d", "--dist", type=float, nargs=1,
		dest="dist", required=False, default=[0.01],
		help="Specify the cutoff distance for clustering in the weighting scheme (default=0.01).")

	# SVM
	parser.add_argument("-c", "--soft_margin", type=float, nargs=1,
		dest="c", required=False, default=[1.0],
		help="The soft margin for the SVM (default=1.0).")
	parser.add_argument("-x", "--xval", type=int, nargs="?",
		dest="xval", required=False, const=5, default=None,
		help="Performs cross validation of training set. Specify folds over 10 repetitions (default=5).")
	parser.add_argument("-e", "--kernel", nargs=2,
		dest="kernel", required=False, default=["linear",0],
		help="Specify kernel to be used and sigma if applicable (i.e. gaussian) (default='linear', 0).")
	parser.add_argument("-s", "--svs", action="store_true",
		dest="svs", required=False,
		help="Show support vectors.")

	return parser

if __name__ == '__main__':
	# Get args
	parser = get_args()
	args = parser.parse_args()
	# Print detail
	mini = args.mini
	### Load training set and make features ###
	gta_file = args.gta[0]
	virus_file = args.virus[0]
	kmer_size = args.kmer[0]
	# Load profiles
	gta_profs = Loader.load(gta_file, "GTA")
	viral_profs = Loader.load(virus_file, "virus")
	# Make features
	feats = Feature(gta_profs.profiles + viral_profs.profiles)
	feats.make_kmer_dict(kmer_size)
	feats.kmer_feat()
	if args.pseaac:
		feats.pseaac(lam=int(args.pseaac), weight=PSE_WEIGHT)

	# Weight if needed
	if args.weight:
		# Get distance threshold
		d = args.dist[0]
		# Get cluster type
		cluster_type = args.cluster_type[0]
		# Weight GTA
		pairwiseGTA = Weight.load(args.weight[0])
		GTA_weight = Weight(gta_profs, pairwiseGTA)
		GTA_clusters = GTA_weight.cluster(cluster_type, d)
		GTA_weight.weight(GTA_clusters)
		# Weight Virus
		pairwiseViral = Weight.load(args.weight[1])
		virus_weight = Weight(viral_profs, pairwiseViral)
		virus_clusters = virus_weight.cluster(cluster_type, d)
		virus_weight.weight(virus_clusters)

	# Create SVM
	c = args.c[0]
	kernel = args.kernel[0]
	sigma = float(args.kernel[1])

	svm = SVM(gta_profs, viral_profs, c, kernel, sigma)

	# Print support vectors
	if args.svs:
		svm.show_svs()

	# Xval
	if args.xval:
		nfolds = args.xval
		result = svm.xval(nfolds, NREPS)
		if mini:
			print("GTA Correct\tViral Correct")
			print("%.2f\t%.2f" % (result[0], result[1]))
		else:
			print("We correctly classified (on average) %.2f/%d GTA and %.2f/%d Viral genes." 
			% (result[0], len(gta_profs), result[1], len(viral_profs)))
	
	else: # Otherwise classify test set
		# Make sure queries set
		if args.queries == None:
			print("The query file was not specified. Please declare queries using -q.")
		else: # All good
			# Load test set
			test_file = args.queries[0]
			test_profs = Loader.load(test_file)
			# Make features
			feats.kmer_feat(test_profs)
			# Classify
			svm.predict(test_profs)
			# Print results
			if mini:
				for profile in test_profs:
					print("Gene\t\tClass")
					print(">%s\t%s" % (55, profile.name, profile.label))
			else:
				print("%-*s%-*s%-*s" % (55, "Gene", 10, "Score", 5, "Classification"))
				for profile in test_profs:
					print(">%-*s%-*f%-*s" % (55, profile.org_name, 10, profile.score, 5, profile.label))