"""
File name: Loader.py
Date created: 5/1/2016
Date last modified: 5/1/2016
Python version: 2.7
Description: A class for generating
	and saving features for a given
	set of profiles.
"""

###############
### IMPORTS ###
###############
from Profile import Profile
from Loader import Loader
import os
import numpy as np

############
### CODE ###
############
class Feature:
	"""Contains all information pertaining to a given gene."""

	def __init__(self, profiles):
		""" Init for Feature. Takes Loader of the
			profiles.
		"""
		self.profiles = profiles
		self.kmer_dict = {}

	@staticmethod
	def load(filename):
		"""Loads dictionary for kmers
		Args:
		  file (string): the name of the file to be loaded
		Returns:
			A Feature Object
		Example:
		"""
		pass

	def make_kmer_dict(k):
		"""Generates dictionary of kmers
			from profile sequences.
		Args:
			k (int): length of kmers
		Returns:
			Bool: Success of function.
		"""
		# Make sure protein sequences loaded
		if profiles[0].prot_seq == None:
			print("Cannot generate kmer features. No protein sequences loaded into Profiles.")
			return False
		# Loop through profiles
		for profile in self.profiles:
			# Loop through profile sequence
			for i in range(len(profile.prot_seq) - k + 1):
				kmer = profile.prot_seq[i:i+k] 
				# Try adding kmer to dictionary
				if kmer not in self.kmer_dict:
					self.kmer_dict[kmer] = 0
				else:
					self.kmer_dict[kmer] += 1

		return True

	def kmer_feat():
		"""Makes feature set using the 
			kmer_dict and appends features
			in respective profiles.
		Args: None
		Returns:
			Bool: Success of function.
		"""
		# Make sure kmers in kmer_dict
		if len(self.kmer_dict) == 0:
			print("No kmers in dictionary. Call make_kmer_dict(k) first.")
			return False

		# Loop through profiles
		for profile in self.profiles:
			# Make temp dictionary?
			temp_kmer_dict = clone(self.kmer_dict)
			# Run through sequence of profile and get kmer count
			for i in range(len(profile.prot_seq) - k + 1):
				kmer = profile.prot_seq[i:i+k] 
				# Try adding kmer to dictionary
				if kmer not in self.kmer_dict:
					self.kmer_dict[kmer] = 0


############################################
############ DANS CODE TO USE ##############
############################################
import Sequence
import numpy as np
from scipy.stats import norm
from collections import defaultdict

"""
Class for converting character sequences into feature vectors using the Spectrum Kernel approach
"""

class Vectorizer:

	def __init__(self, k, var=True):
		self.k = k 		# length of k-mers 	
		self.var = var  # boolean indicating whether to vary size of k-mers from 1 to k or fix to size at k.
		self.bag = None	# set of k-mers comprising "spectrum"

	"""
	Scan character sequence with window of size k, counting the occurences of substrings of length k.

	Input:
	seq = character string

	Output:
	bag = defaultdict mapping k-mers to their frequency within seq
	"""

	def get_kmer_counts(self, seq):
		bag = defaultdict(int)
		if self.var == True:
			K = range(1, self.k+1)
		else:
			K = [k]
		for k in K:
			for x in range(len(seq)-k+1):
				kmer = seq[x:x+k]
				bag[kmer] += 1
		return bag

	"""
	Converts character sequence into feature vector.

	Input:
	seq = character string

	Output:
	feat_vec = list containing k-mer frequency counts for input sequence
	corres_kmers = list containing k-mer strings corresponding to frequencies in feat_vec
	"""

	def vectorize(self, seq):
		counts = self.get_kmer_counts(seq) # get dictionary of k-mer counts
		feat_vec = []
		corres_kmers = []

		# create feature vector by referring to k-mer dictionary
		for word in self.bag:
			feat_vec.append(counts[word])
			corres_kmers.append(word)
		return feat_vec, corres_kmers

	"""
	Identifies "sparse" features and returns a list of the indicies of these features

	Input:
	train_set = 2-D numpy array of feature vectors in training set
	thr = threshold for determining sparsity. If thr=4, any k-mer appearing less than 4 times in the entire
		  training set will be considered a "sparse" feature

	Output:
	to_prune = list of sparse feature indicies 
	"""

	def get_sparse_features(self, train_set, thr):
		zipped = zip(*train_set) 
		to_prune = [i for i in range(len(train_set[0])) if sum(zipped[i]) < thr]		
		to_prune.sort(reverse=True) # sort sparse feature indicies from greatest to least so that you can delete features
									# in order, without affecting the indicies of other sparse features
		return to_prune

	"""
	Returns indicies of n best features according to Bi-Normal Selection (BNS) metric.

	Input:
	X = 2-D Numpy array of feature vectors for training sequences
	y = 1-D Numpy array of class labels for training sequences
	n = number of features to select using BNS

	Output:
	n_best_feats_idx = indicies of the n best features
	"""

	def select_features(self, X, y, n=1000):		
		zipped = zip(*X)
		scores = []
		pos_scores = []
		neg_scores = []
		
		# Evaluate each feature
		for col in zipped:
			col = np.array(list(col))
			
			thr = np.average(col)	# set threshold value for feature 
			
			# get positive and negative rates
			pos = col[y==1]
			tp = pos[pos>thr]				
			neg = col[y!=-1]
			fp = neg[neg>thr]
			prev_pos = float(len(tp))/len(pos)
			prev_neg = float(len(fp))/len(neg)

			# compute BNS score
			bns_score = abs(norm.ppf(prev_pos) - norm.ppf(prev_neg))
			scores.append(bns_score)

			if prev_pos > prev_neg:
				pos_scores.append(bns_score)
				neg_scores.append('nope')
			else:
				pos_scores.append(bns_score)
				neg_scores.append('nope')

		# get the indicies of the n best features, and sort them from most best to least best
		n_best_feats_idx = np.array(scores).argsort()[-n:]
		n_best_feats_idx = n_best_feats_idx[::-1]

		# extract the indicies corresponding to GTA-indicative features and Virus-indicative features
		pos_feats_idx = [i for i in n_best_feats_idx if type(pos_scores[i]) is not str]
		neg_feats_idx = [i for i in n_best_feats_idx if type(neg_scores[i]) is not str]
		
		return n_best_feats_idx








