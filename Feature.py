"""
File name: Loader.py
Date created: 01/05/2016
Date last modified: 04/07/2016
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

	def kmer_select(self, n):
		""" Selects the top n most
			informative features from the 
			kmer set. Called after kmers
			are generated... Gonna be complex.
		"""
		keepers = []

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

# Command-line driver
if __name__ == '__main__':
	def standard_convert(param_dict):
		""" Standard conversion for given AA calculated from
			the sum of the 20 amino acids
		"""
		dicter = {}
		for aa in AA_LETTERS:
			H = param_dict[aa]
			Ha = sum([param_dict[key] for key in AA_LETTERS])/20.0
			numerator = H - Ha
			denominator = math.sqrt(sum([(param_dict[key] - Ha)**2 for key in AA_LETTERS])/20.0)
			dicter[aa] = numerator / denominator
		return dicter

	print("HYDROPHIL_STANDARD = ",standard_convert(HYDROPHIL_SCORES))
	print("HYDROPHOB_STANDARD = ",standard_convert(HYDROPHOB_SCORES))
	print("SIDECHAIN_STANDARD = ",standard_convert(SIDECHAIN_MASS))
	# # Define arg parser
	# parser = argparse.ArgumentParser(description="Gene featuring.")

	# ### Define Args ###

	# # Main
	# parser.add_argument("-g", "--GTA", type=str, nargs=1,
	# 	dest="gta", required=True,
	# 	help="The .faa or .fna training file for GTA genes.")
	# parser.add_argument("-v", "--virus", type=str, nargs=1,
	# 	dest="virus", required=True,
	# 	help="The .faa or .fna training file for viral genes.")
	# parser.add_argument("-k", "--kmer", type=int, nargs=1,
	# 	dest="kmer", required=False, default=[4],
	# 	help="The kmer size needed for feature generation (default=4).")

	# args = parser.parse_args()
	
	# ### Load training set and make features ###
	# gta_file = args.gta[0]
	# virus_file = args.virus[0]
	# kmer_size = args.kmer[0]
	# # Load profiles
	# gta_profs = Loader.load(gta_file, "GTA")
	# viral_profs = Loader.load(virus_file, "virus")
	# # Make features
	# feats = Feature(gta_profs.profiles + viral_profs.profiles)
	# feats.pseaac(lam=3, weight=0.05)
	# feats.make_kmer_dict(kmer_size)
	# feats.kmer_feat(gta_profs.profiles + viral_profs.profiles)

	# ## Test Code ###
	# # misclassified viral trainer: 77864702
	# # correct viral trainer: 712912299, 168495136, 23505451
	# # misclassified gta trainers: 85375630, 222107096
	# # correct gta trainer: 126461901

	# featureSet = np.array([x.features for x in gta_profs] + [y.features for y in viral_profs])
	# featureSum = np.sum(featureSet, 0)

	# v1 = viral_profs.profileDict["77864702"].features
	# commonv1 = [1 if featureSum[i] == v1[i] else 0 for i in range(len(v1))]
	# print("%s has %d/%d unique kmers in the kmer bag." % (viral_profs.profileDict["77864702"].org_name, sum(commonv1), len(commonv1)))

	# v2 = viral_profs.profileDict["23505451"].features
	# commonv2 = [1 if featureSum[i] == v2[i] else 0 for i in range(len(v2))]
	# print("%s has %d/%d unique kmers in the kmer bag." % (viral_profs.profileDict["23505451"].org_name, sum(commonv2), len(commonv2)))

	# v3 = viral_profs.profileDict["168495136"].features
	# commonv3 = [1 if featureSum[i] == v3[i] else 0 for i in range(len(v3))]
	# print("%s has %d/%d unique kmers in the kmer bag." % (viral_profs.profileDict["168495136"].org_name, sum(commonv3), len(commonv3)))

	# g1 = gta_profs.profileDict["85375630"].features
	# commong1 = [1 if featureSum[i] == g1[i] else 0 for i in range(len(g1))]
	# print("%s has %d/%d unique kmers in the kmer bag." % (gta_profs.profileDict["85375630"].org_name,sum(commong1), len(commong1)))

	# g2 = gta_profs.profileDict["222107096"].features
	# commong2 = [1 if featureSum[i] == g2[i] else 0 for i in range(len(g2))]
	# print("%s has %d/%d unique kmers in the kmer bag." % (gta_profs.profileDict["222107096"].org_name,sum(commong2), len(commong2)))

	# g3 = gta_profs.profileDict["126461901"].features
	# commong3 = [1 if featureSum[i] == g3[i] else 0 for i in range(len(g3))]
	# print("%s has %d/%d unique kmers in the kmer bag." % (gta_profs.profileDict["126461901"].org_name,sum(commong3), len(commong3)))
	# # azosprill and burkholderia commonality
	# commons = [1 if  v1[i] > 0 and v3[i] > 0 else 0 for i in range(len(v1))]
	# print("%s and %s have %d/%d kmers in common." % 
	# 	(viral_profs.profileDict["77864702"].org_name,
	# 	viral_profs.profileDict["168495136"].org_name, sum(commons), len(commons)))
	# # azosprill and lactobac commonality
	# commonss = [1 if  v2[i] > 0 and v3[i] > 0 else 0 for i in range(len(v1))]
	# print("%s and %s have %d/%d kmers in common." % 
	# 	(viral_profs.profileDict["23505451"].org_name,
	# 	viral_profs.profileDict["168495136"].org_name, sum(commonss), len(commonss)))
	# # burkholderia and lactobac commonality
	# commonsss = [1 if  v1[i] > 0 and v2[i] > 0 else 0 for i in range(len(v1))]
	# print("%s and %s have %d/%d kmers in common." % 
	# 	(viral_profs.profileDict["77864702"].org_name,
	# 	viral_profs.profileDict["23505451"].org_name, sum(commonsss), len(commonsss)))
	# # burkholderia and lactobac commonality
	# commonssss = [1 if  g1[i] > 0 and g2[i] > 0 else 0 for i in range(len(v1))]
	# print("%s and %s have %d/%d kmers in common." % 
	# 	(gta_profs.profileDict["85375630"].org_name,
	# 	gta_profs.profileDict["222107096"].org_name, sum(commonssss), len(commonssss)))

	# print(np.count_nonzero(gta_profs.profileDict["222107096"].features))
