"""
File name: Profile.py
Date created: 22/10/2015
Date last modified: 22/10/2015
Python version: 2.7
Description: A data structure for
	gene data. Used by Loader in 
	GTA_Tool suite.
"""

############
### CODE ###
############
class Profile:
	"""Contains all information pertaining to a given gene."""

	def __init__(self, name, label=None, dna_seq=None, prot_seq=None, weight=None, features=None, score=None):
		"""A Profile with at least one initialized variable
		Args:
		  name (string): name of the gene
		  label (string): label of gene (optional)
		  dna_seq (string): DNA sequence from FASTA (optional)
		  prot_seq (string): AA seqeuence from FASTA (optional)
		  weight (float): given weight when training with profile (optional)
		  features: a sequence of features derived from Feature (optional)
		"""
		self.name = name
		self.label = label
		self.dna_seq = dna_seq
		self.prot_seq = prot_seq
		self.weight = weight
		self.features = features
		self.score = score

	def __str__(self):
		"""The name and class of the Profile"""
		out = self.name
		if self.label:
			out = out + ": " + self.label
		if self.score:
		  out = out + "\t" + score
		return out

	def __getitem__(self, pos):
		"""The feature at given postion
		Returns:
		  Feature (float)
		"""
		return self.features[pos]

	def __len__(self):
		"""The length of features
		Returns:
		  int
		"""
		return len(self.features)
