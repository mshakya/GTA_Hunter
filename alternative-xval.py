###############
### IMPORTS ###
###############
from Loader import Loader
from Weight import Weight
from Feature import Feature
import argparse
import sys
import random
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import MultinomialNB

#################
### CONSTANTS ###
#################
MINI = False
NFOLDS = 10 
NREPS = 10 
PSE_WEIGHT = 0.05 # for pseaac feature
# Specifics
GENE = "15"
K = 4
LAM = 15
# Dir
DIR = "../data/2016-3-28/"
OUT = "../results/2016-5-5/xval.txt"
GTA_PATH = DIR + "gta/" + GENE + "_gta.faa"
VIRAL_PATH = DIR + "viral/" + GENE + "_viral.faa"


############
### CODE ###
############
def xval(predictor, train0, train1, nfold=5, nrep=10):
	"""n-fold cross validation of
		the test set. 
		Input:
			n (int): number of folds for xval
		Returns:
			(fpr, fnr): false positive and false 
				negative rates from the n xvals
	"""
	# Set class labels
	label0 = train0.label
	label1 = train1.label
	profiles = train0.profiles + train1.profiles
	# keep track of label classification
	score0 = 0.0
	score1 = 0.0
	gta_as_phage = []
	phage_as_gta = []

	# repeat xval results nrep times
	for i in range(nrep):
		if not MINI:
			sys.stdout.flush()
			sys.stdout.write("Starting rep: %d\r" % (i+1))
		# randomly sort profiles
		random.shuffle(profiles)
		# split into folds
		split = [profiles[i::nfold] for i in range(nfold)]
		# cross val
		for j in range(nfold):
			# Build train and test sets
			train_fold = np.array([x for sublist in (split[:j]+split[j+1:]) for x in sublist])
			test_fold = split[j]
			trainX = np.array([x.features for x in train_fold])
			testX = np.array([x.features for x in test_fold])
			trainY = np.array([-1.0 if y.label == label0 else 1.0 for y in train_fold])
			testY = np.array([-1.0 if y.label == label0 else 1.0 for y in test_fold])
			testNames = np.array([x.org_name for x in test_fold])
			# Get training set weights

			# Evaluate results
			predictor.fit(trainX, trainY)
			for r in range(len(testX)):
				# Positive product is correct classification
				if predictor.predict([testX[r]]) * testY[r] > 0:
					# Update label0 if negative, label1 otherwise
					if testY[r] < 0:
						score0 += 1
					else:
						score1 += 1
				else: # predicted incorrectly
					if testY[r] > 0: #virus as GTA
						phage_as_gta.append(testNames[r])
					else: #gta as virus
						gta_as_phage.append(testNames[r])

	# if not MINI:
	# 	print("\nPhages (%d) misclassified over %d reps: %s" % (len(phage_as_gta), nrep, phage_as_gta))
	# 	print("\nGTA (%d) misclassified over %d reps: %s\n" % (len(gta_as_phage), nrep, gta_as_phage))

	return (score0/nrep, score1/nrep)

if __name__ == '__main__':
	# Load profiles
	gta_profs = Loader.load(GTA_PATH, "GTA")
	viral_profs = Loader.load(VIRAL_PATH, "virus")
	# Make features
	feats = Feature(gta_profs.profiles + viral_profs.profiles)
	# kmer
	feats.make_kmer_dict(K)
	feats.kmer_feat()
	# pseaac
	feats.pseaac(lam=LAM, weight=PSE_WEIGHT)
	# physicochem
	feats.physicochem()

	# Xval
	# predictor = KNeighborsClassifier(n_neighbors=10)
	predictor = MultinomialNB()
	result = xval(predictor, gta_profs, viral_profs, NFOLDS, NREPS)
	if MINI:
		print("GTA Correct\tViral Correct")
		print("%.2f\t%.2f" % (result[0], result[1]))
	else:
		print("We correctly classified (on average) %.2f/%d GTA and %.2f/%d Viral genes." 
		% (result[0], len(gta_profs), result[1], len(viral_profs)))