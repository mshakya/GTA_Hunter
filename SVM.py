"""
File name: SVM.py
Date created: 4/12/2015
Date last modified: 5/1/2016
Python version: 2.7
Description: Support Vector Machine
	class that can train, test, and
	cross validate on an array of Profiles.
"""

from Profile import Profile
from Loader import Loader
import cvxopt
import numpy as np
import sys

tau = 1e-5

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

###########
### SVM ###
###########
class SVM:
	"""Support Vector Machine Class.
		Contains train and predict inner
		classes.
	"""
	#############
	### TRAIN ###
	#############
	class SVMTrain:
		"""Trainer for SVM. Takes kernel and c
			(soft margin) and allows for training
			on various training sets.
		"""

		def __init__(self, kernel, c, weights):
			"""Initializes trainer
				Input:
					kernel (lambda function): the kernel to be used
					c (float): soft margin var
				Returns:
					None
			"""
			self.kernel = kernel
			self.c = c
			self.weights = weights

		def train(self, trainX, trainY):
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
			lagrange = self.compute_lagrange(trainX, trainY)

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

		def compute_lagrange(self, trainX, trainY):
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

			# Soft Margin
			elif self.weights == None:
				G = cvxopt.matrix(np.vstack((-np.eye(m_profiles), np.eye(m_profiles))))
				h = cvxopt.matrix(np.hstack((np.zeros(m_profiles), self.c * np.ones(m_profiles))))

			# Weighted Soft Margin
			else:
				G = cvxopt.matrix(np.vstack((-np.eye(m_profiles), np.eye(m_profiles))))
				h = cvxopt.matrix(np.hstack((np.zeros(m_profiles), self.c * self.weights)))

			##############

			# Other Pagams
			P = cvxopt.matrix(np.outer(trainY, trainY) * K)
			q = cvxopt.matrix(-1 * np.ones(m_profiles))
			A = cvxopt.matrix(trainY, (1,m_profiles))
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
			# refine this
			result = self._bias
			for z_i, x_i, y_i in zip(self.alphas,
								 self.svs,
								 self.svsLabels):
				result += z_i * y_i * self.kernel(x_i, x)

			return np.sign(result).item()

	def __init__(self, train0, train1, c, kernel, sigma=None, xval=False):
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
		self.trainY = np.array([0] * len(train0) + [1] * len(train1))
		# Get training set weights
		if train0[0].weight != None and train1[0].weight != None:
			self.weights = np.array([x.weight for x in train0] + [y.weight for y in train1])
		else:
			self.weights = None
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
		# Build Predictor if not cross validating training set
		if not xval:
			self.predictor = self.SVMTrain(self.kernel, self.c, self.weights).train(self.trainX, self.trainY)


	def __init__(self, *args):
		""" Overloaded method.
			Do not call directly, use the load 
			method instead. Only inits predictor.
			Example:
				SVM = SVM.load("file.svm")
		"""
		# Deconstruct args
		kernel = args[0]
		bias = args[1]
		alphas = args[2]
		svs = args[3]
		svsLabels = args[4]
		# Build Predictor
		self.predictor = self.SVMTrain(kernel, bias, alphas, svs, svsLabels)
		

	@staticmethod
	def load(filename):
		"""Loads given the trained params into the SVM
		Args:
		  file (string): the name of the file to be loaded
		Returns:
			SVM object with predictor
		Example:
			SVM = SVM.load("file.svm")
		"""
		# Code goes here.
		pass
		# Return SVM model
		args = (kernel, bias, alphas, svs, svsLabels)
		return SVM(args)

	def xval(self, n):
		"""n-fold cross validation of
			the test set. 
			Input:
				n (int): number of folds for xval
			Returns:
				(fpr, fnr): false positive and false 
					negative rates from the n xvals
		"""

# Command-line driver -- just some hard-coded test cases -- add your own if you want
if __name__ == '__main__':
	# Get classify or xval
	command = argv[1]
	profiles = ExpressionProfile.read(argv[2])
	if command=='classify':
		train0 = input()
	elif command=='xval':
		print("xvaling")
	else:
		print("Unknown command. Accepts the following: classify, xval")
