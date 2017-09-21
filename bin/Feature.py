"""
File name: Loader.py
Date created: 01/05/2016
Date last modified: 05/31/2016
Python version: 3.5.1
Description: A class for generating
	and saving features for a given
	set of profiles.
"""

###############
### IMPORTS ###
###############
from Profile import Profile
from Loader import Loader
import argparse
import os
import numpy as np
import math

#################
### CONSTANTS ###
#################
AA_LETTERS = ['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 
				'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T']

# Scores and masses taken from http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/ParaValue.htm table
HYDROPHOB_SCORES = {'A': 0.62, 'C': 0.29, 'D': -0.9, 'E': -0.74, 'F': 1.19, 
					'G': 0.48, 'H': -0.4, 'I': 1.38, 'K': -1.5, 'L': 1.06, 
					'M': 0.64, 'N': -0.78, 'P': 0.12, 'Q': -0.85, 'R': -2.53, 
					'S': -0.18, 'T': -0.05, 'V': 1.08, 'W': 0.81, 'Y': 0.26}

HYDROPHIL_SCORES = {'A': -0.5, 'C': -1.0, 'D': 3.0, 'E': 3.0, 'F': -2.5, 
					'G': 0.0, 'H': -0.5, 'I': -1.8, 'K': 3.0, 'L': -1.8, 
					'M': -1.3, 'N': 0.2, 'P': 0.0, 'Q': 0.2, 'R': 3.0, 
					'S': 0.3, 'T': -0.4, 'V': -1.5, 'W': -3.4, 'Y': -2.3}

SIDECHAIN_MASS = {'A': 15.0, 'C': 47.0, 'D': 59.0, 'E': 73.0, 'F': 91.0, 
					'G': 1.0, 'H': 82.0, 'I': 57.0, 'K': 73.0, 'L': 57.0, 
					'M': 75.0, 'N': 58.0, 'P': 42.0, 'Q': 72.0, 'R': 101.0, 
					'S': 31.0, 'T': 45.0, 'V': 43.0, 'W': 130.0, 'Y': 107.0}

# Standard conversion of scores and masses above, according to http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/type1.htm
HYDROPHIL_STANDARD =  {'A': -0.15187800657861578, 'Y': -1.1111075218119784, 'Q': 0.2211556937899142, 'L': -0.8446548786915999, 'N': 0.2211556937899142, 'C': -0.4183306496989943, 'H': -0.15187800657861578, 'R': 1.7132904952640338, 'I': -0.8446548786915999, 'D': 1.7132904952640338, 'S': 0.27444622241398986, 'T': -0.09858747795454006, 'G': 0.11457463654176277, 'V': -0.6847832928193728, 'W': -1.6973033366768113, 'F': -1.21768857906013, 'E': 1.7132904952640338, 'P': 0.11457463654176277, 'K': 1.7132904952640338, 'M': -0.5782022355712214}

HYDROPHOB_STANDARD =  {'A': 0.6362505881446506, 'Y': 0.26681476277033744, 'Q': -0.8722790321337952, 'L': 1.0877832636021447, 'N': -0.8004442883110121, 'C': 0.2976010815515302, 'H': -0.41048425041590364, 'R': -2.5963128838805902, 'I': 1.4161706639348675, 'D': -0.9235895634357832, 'S': -0.18471791268715662, 'T': -0.051310531301987934, 'G': 0.49258110049908443, 'V': 1.10830747612294, 'W': 0.831230607092205, 'F': 1.2211906449873133, 'E': -0.7593958632694218, 'P': 0.12314527512477112, 'K': -1.5393159390596387, 'M': 0.6567748006654459}

SIDECHAIN_STANDARD =  {'A': -1.5919364305641373, 'Y': 1.4624567208832167, 'Q': 0.3004593263108537, 'L': -0.19753955707730178, 'N': -0.16433963151809142, 'C': -0.5295388126694055, 'H': 0.6324585819029575, 'R': 1.2632571675279545, 'I': -0.19753955707730178, 'D': -0.13113970595888105, 'S': -1.0607376216167714, 'T': -0.5959386637878262, 'G': -2.0567353883930823, 'V': -0.662338514906247, 'W': 2.226055008745055, 'F': 0.9312579119358507, 'E': 0.3336592518700641, 'P': -0.6955384404654573, 'K': 0.3336592518700641, 'M': 0.40005910298848485}

# Physicochemical properties of amino acids from Kaundal et al. Bioinformatics 2013
PHYSICOCHEM = {'A': [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
				'C': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1],
				'D': [1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
				'E': [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
				'F': [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
				'G': [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
				'H': [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
				'I': [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
				'K': [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
				'L': [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
				'M': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
				'N': [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0],
				'P': [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
				'Q': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0],
				'R': [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
				'S': [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
				'T': [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
				'V': [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
				'W': [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
				'Y': [0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0]}

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

	def make_kmer_dict(self, k, n=0):
		"""Generates dictionary of kmers
			from profile sequences.
		Args:
			k (int): length of kmers
		Returns:
			Bool: Success of function.
		"""
		# Make sure protein sequences loaded
		if self.profiles[0].prot_seq == None:
			print("Cannot generate kmer features. No protein sequences loaded into Profiles.")
			return False
		# Loop through profiles
		for profile in self.profiles:
			# Loop through profile sequence
			profile_kmers = [0] * (len(profile.prot_seq) - k + 1)
			for i in range(len(profile.prot_seq) - k + 1):
				kmer = profile.prot_seq[i:i+k] 
				# Try adding kmer to dictionary
				if kmer not in self.kmer_dict:
					self.kmer_dict[kmer] = 1
				elif kmer not in profile_kmers:
					self.kmer_dict[kmer] += 1
				profile_kmers[i] = kmer

		# Remove unique features (kmer count == 1)
		self.kmer_dict = {kmer: value for kmer, value in self.kmer_dict.items() if value > 1}


		# save k
		self.k = k
		return True

	def kmer_feat(self, profiles=None):
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
		# Set profiles
		if profiles == None:
			profiles = self.profiles
		# Loop through profiles
		for profile in profiles:
			# Make temp dictionary inititing all keys at 0
			profile_kmer_dict = dict.fromkeys(self.kmer_dict, 0)
			# Run through sequence of profile and get kmer count
			for i in range(len(profile.prot_seq) - self.k + 1):
				kmer = profile.prot_seq[i:i+self.k] 
				# Try adding kmer to dictionary
				if kmer in profile_kmer_dict:
					profile_kmer_dict[kmer] += 1
			# Rip dictionary  and add features (values) to profile
			feats = [profile_kmer_dict[key] for key in self.kmer_dict]
			if profile.features:
				profile.features += feats 
			else:
				profile.features = feats

	def aa_score(self, aa1, aa2):
		""" Returns score of two amino acids for
			use in pseaac calculations. Incorporates
			the following scoring parameters:
				hydrophobicity value
				hydrophilicity value
				side chain mass
		"""
		hydrophobicity = (HYDROPHOB_STANDARD[aa2] - HYDROPHOB_STANDARD[aa1])**2
		hydrophilicity = (HYDROPHIL_STANDARD[aa2] - HYDROPHIL_STANDARD[aa1])**2
		side_chain = (SIDECHAIN_STANDARD[aa2] - SIDECHAIN_STANDARD[aa1])**2

		return (hydrophilicity + hydrophobicity + side_chain) / 3

	def pseaac(self, lam=3, weight=0.05, profiles=None):
		""" Pseudo amino acid composition.
			Appends PseAA composition features
			to profile feature set. Unlike traditional
			AAC, PseAAC preserves some sequence-order
			information.
		"""
		# Set profiles
		if profiles == None:
			profiles = self.profiles
		# Loop through profiles
		for profile in profiles:
			# Calculate the first 20 features (AAs)
			p20 = {aa: 0.0 for aa in AA_LETTERS}
			for aa in profile.prot_seq:
				p20[aa] += 1
			# Turn p20 into numeric list, keeping order constant
			p20 = [value for (key, value) in sorted(p20.items())]
			# Normalize
			p20 = [aa/sum(p20) for aa in p20]
			# Calculate the tier scores (lambda)
			p20plus = [0.0] * lam
			for k in range(1, lam+1):
				J = 0
				# Sum kth most contiguous residues
				for i in range(len(profile.prot_seq) - k):
					J += self.aa_score(profile.prot_seq[i], profile.prot_seq[i+k])
				# Average
				J = J / (len(profile.prot_seq) - k)
				# Weight and store
				p20plus[k-1] = J * weight
			# Normalize
			normalizer = sum(p20) + sum(p20plus)
			# Create feature set
			feats = [x / normalizer * 100 for x in p20 + p20plus]
			# Store
			if profile.features:
				profile.features += feats
			else:
				profile.features = feats

	def physicochem(self, profiles=None):
		""" Physicochemical properties of amino
			acids. Each protein generates a 19
			dimensional feature vector consiting
			of normalized frequencies of each property.
		"""
		# Set profiles
		if profiles == None:
			profiles = self.profiles
		# Loop through profiles
		for profile in profiles:
			# Initialize features
			feats = [0.0] * 19
			# Iterate over protein sequence
			for aa in profile.prot_seq:
				feats = np.add(feats, PHYSICOCHEM[aa])
			# Normalize by sequence length
			feats = [x / len(profile.prot_seq) for x in feats]
			# Store
			if profile.features:
				profile.features += feats
			else:
				profile.features = feats
