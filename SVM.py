"""
File name: SVM.py
Date created: 4/12/2015
Date last modified: 5/1/2016
Python version: 3.5
Description: Support Vector Machine
	class that can train, test, and
	cross validate on an array of Profiles.
"""

###############
### IMPORTS ###
###############
from Profile import Profile
from Loader import Loader
import numpy as np
import cvxopt
import argparse
import sys
import random

# SV Cutoff and CVXOPT print control
tau = 1e-5
cvxopt.solvers.options['show_progress'] = False
mini = True

##############
### KERNEL ###
##############
class Kernel:
	"""Various kernels to be used by SVM.
		Only contains staticmethods
		Example:
			x = Kernel.linear(x,y)
	"""

	# NOTE: all methods use return lambda, to allow 
	# for simple x,y function handling in trainer and
	# predictor methods.

	@staticmethod
	def linear():
		return lambda x, y: np.inner(x,y)

	@staticmethod
	def gaussian(sigma):
		return lambda x, y: np.exp(-np.sqrt(np.linalg.norm(x-y)**2 / (2*sigma**2)))
#############
### TRAIN ###
#############
class SVMTrain:
	"""Trainer for SVM. Takes kernel and c
		(soft margin) and allows for training
		on various training sets.
	"""

	def __init__(self, kernel, c):
		"""Initializes trainer
			Input:
				kernel (lambda function): the kernel to be used
				c (float): soft margin var
			Returns:
				None
		"""
		self.kernel = kernel
		self.c = c

	def train(self, trainX, trainY, weights):
		"""Trains SVM model for given training
			set. Calls functions to compute params.
			Input: 
				trainX[n,m] (object): a numpy array of Profiles
					and features
				trainY[1,n] (object): a numpy vector assigning each
					Profile in trainX to a class
			Returns:
				SVMPredict: an initialized predictor
		"""
		# get lagrange multipliers
		lagrange = self.compute_lagrange(trainX, trainY, weights)

		# find support vectors and alphas
		svinds = lagrange > tau
		alphas = lagrange[svinds]
		svs = trainX[svinds]
		svsLabels = trainY[svinds]

		# compute bias by calculating error from bias=0
		bias = np.mean([y_i - SVMPredict(self.kernel, 0.0, alphas, svs, 
			svsLabels).predict(x_i) for (y_i, x_i) in zip(svsLabels, svs)])
		
		return SVMPredict(self.kernel, bias, alphas, svs, svsLabels)

	def compute_k(self, trainX):
		"""Computes K, the gram matrix X'X, using
			the kernel.
			Input:
				trainX[m,n]: numpy matrix of profiles and features
			Returns:
				K[m,m]: numpy matrix for use in quadprog
		"""
		# Get dimensions of X
		# Row = a profile (m)
		# Column = a feature (n)
		m_profiles, n_features = trainX.shape # returns row, column tuple

		# Init K
		K = np.zeros((m_profiles,m_profiles))

		# Compute each i,j in K
		for i in range(m_profiles):
			for j in range(m_profiles):
				K[i,j] = self.kernel(trainX[i], trainX[j])

		return K

	def compute_lagrange(self, trainX, trainY, weights):
		"""Computes lagrange multipliers by computing
			various variables needed for quadprog.
			Input:
				trainX[m,n] (object): a numpy array of Profiles
					and features
				trainY[1,m] (object): a numpy vector assigning each
					Profile in trainX to a class
			Returns:
				lagrange multipliers GIVE MORE DETAIL ABOUT WHAT THIS IS 
		"""
		# Get dimensions of X
		# Row = a profile (m)
		# Column = a feature (n)
		m_profiles, n_features = trainX.shape # returns row, column tuple

		# Compute K
		K = self.compute_k(trainX)

		##############
		### Margin ###
		##############

		# Hard Margin
		if not self.c:
			G = cvxopt.matrix(np.diag(np.ones(m_profiles) * -1))
			h = cvxopt.matrix(np.zeros(m_profiles))

		# Weighted Soft Margin ( if no weights, weights = ones )
		else:
			G = cvxopt.matrix(np.vstack((-np.eye(m_profiles), np.eye(m_profiles))))
			h = cvxopt.matrix(np.hstack((np.zeros(m_profiles), self.c * weights)))

		##############

		# Other Pagams
		P = cvxopt.matrix(np.outer(trainY, trainY) * K)
		q = cvxopt.matrix(-1 * np.ones(m_profiles))
		A = cvxopt.matrix(trainY, (1, m_profiles))
		b = cvxopt.matrix(0.0)
		# get z (lagrange multipliers)
		z = cvxopt.solvers.qp(P, q, G, h, A, b)

		return np.ravel(z['x'])

###############
### PREDICT ###
###############
class SVMPredict:
	"""Predictor for SVM, takes params and 
		computes label for test cases.				
	"""

	def __init__(self, kernel, bias, alphas, svs, svsLabels):
		"""Initializes predictor
			Input:
				Kernel (staticmethod): The kernel to be used
				Alphas (numpy vector): support vector lagrange multipliers
				Bias (float): Offest from origin
				Classes (list[0,1] of strings): Names of classes 
			Returns:
				None
		"""
		self.kernel = kernel
		self.bias = bias
		self.alphas = alphas
		self.svs = svs
		self.svsLabels = svsLabels

	def predict(self, x):
		"""Predicts classes of test set
			Input:
				Profile: a single profile with features
			Returns:
				None: Updates .label and 
				.score of profile and has them

				Score (float): score for given profile
		"""
		# refine this?
		result = self.bias
		for z_i, x_i, y_i in zip(self.alphas,
							 self.svs,
							 self.svsLabels):
			result += z_i * y_i * self.kernel(x_i, x)

		return result

###########
### SVM ###
###########
class SVM:
	"""Support Vector Machine Class.
		Contains train and predict inner
		classes.
	"""
	
	def __init__(self, train0, train1, c, kernel, sigma=None):
		""" Main init of SVM.
			Input:
				train: training sets two Loaders of Profile arrays for the two classes
				test: test set of Profiles to be classified after training
		"""
		# Set class labels
		self.label0 = train0.label
		self.label1 = train1.label
		# Create training set
		self.trainX = np.array([x.features for x in train0] + [y.features for y in train1])
		self.trainY = np.array([-1.0] * len(train0) + [1.0] * len(train1))
		# Get training set weights
		if train0[0].weight != None and train1[0].weight != None:
			self.weights = np.array([x.weight for x in train0] + [y.weight for y in train1])
		else:
			if not mini:
				print("No weights found, defaulting to ones.")
			self.weights = np.array([1 for x in train0] + [1 for y in train1])
		# Set soft margin
		self.c = c
		# Establish Kernel
		if kernel == "linear":
			self.kernel = Kernel.linear()
		elif kernel == "gaussian":
			try:
				self.kernel = Kernel.gaussian(float(sigma))
			except (ValueError, TypeError):
				raise Exception("The sigma value could not be processed. Try and float or int.")
		else:
			raise Exception("Kernel type could not be detected. Try the following: linear, gaussian, MORE IF ADDED")
		# Build Predictor
		self.predictor = SVMTrain(self.kernel, self.c).train(self.trainX, self.trainY, self.weights)
		# Save train0 and train1 for xval
		self.profiles = train0.profiles + train1.profiles

	def predict(self, test_profiles):
		# Predict each profile in test set
		for profile in test_profiles:
			# Get feats
			feats = np.array(profile.features)
			# Get score
			score = self.predictor.predict(feats)
			# Get label
			if score < 0:
				label = self.label0
			else:
				label = self.label1
			# Assign to profile
			profile.score = score
			profile.label = label

	def xval(self, nfold=5, nrep=10):
		"""n-fold cross validation of
			the test set. 
			Input:
				n (int): number of folds for xval
			Returns:
				(fpr, fnr): false positive and false 
					negative rates from the n xvals
		"""
		# keep track of label classification
		score0 = 0.0
		score1 = 0.0
		gta_as_phage = []
		phage_as_gta = []

		# repeat xval results nrep times
		for i in range(nrep):
			if not mini:
				sys.stdout.flush()
				sys.stdout.write("Starting rep: %d\r" % (i+1))
			# randomly sort profiles
			random.shuffle(self.profiles)
			# split into folds
			split = [self.profiles[i::nfold] for i in range(nfold)]
			# cross val
			for j in range(nfold):
				# Build train and test sets
				train_fold = np.array([x for sublist in (split[:j]+split[j+1:]) for x in sublist])
				test_fold = split[j]
				trainX = np.array([x.features for x in train_fold])
				testX = np.array([x.features for x in test_fold])
				trainY = np.array([-1.0 if y.label == self.label0 else 1.0 for y in train_fold])
				testY = np.array([-1.0 if y.label == self.label0 else 1.0 for y in test_fold])
				testNames = np.array([x.org_name for x in test_fold])
				# Get training set weights
				if np.any(self.weights != 1):
					weights = np.array([x.weight for x in train_fold])
				else:
					weights = np.array([1 for x in train_fold])
				# evaluate results
				predictor = SVMTrain(self.kernel, self.c).train(trainX, trainY, weights)
				for r in range(len(testX)):
					# Positive product is correct classification
					if predictor.predict(testX[r]) * testY[r] > 0:
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

		if not mini:
			print("\nPhages (%d) misclassified over %d reps: %s" % (len(phage_as_gta), nrep, phage_as_gta))
			print("\nGTA (%d) misclassified over %d reps: %s\n" % (len(gta_as_phage), nrep, gta_as_phage))

		return (score0/nrep, score1/nrep)

	def show_svs(self):
		# Get support vectors
		svs = self.predictor.svs
		# Look up index of support vectors in training set
		svsi = [np.where((self.trainX == sv).all(axis=1))[0].tolist()[0] for sv in svs]
		# Map indices back to profiles
		svsp = [self.profiles[i] for i in svsi]
		# Print results
		print("\nThere are %d support vectors." % len(svs))
		print("The support vectors are as follows:")
		for profile in svsp:
			print(profile)
		print()