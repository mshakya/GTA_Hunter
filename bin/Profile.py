"""
File name: Profile.py
Date created: 10/22/2015
Date last modified: 04/07/2016
Python version: 3.5.1
Description: A data structure for
	gene data. Used by Loader in 
	GTA_Tool suite.
"""

############
### CODE ###
############
class Profile:
	"""Contains all information pertaining to a given gene."""

	def __init__(self, name, org_name, label=None, dna_seq=None, prot_seq=None, weight=None, features=None, score=None):
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
		self.org_name = org_name
		self.label = label
		self.dna_seq = dna_seq
		self.prot_seq = prot_seq
		self.weight = weight
		self.features = features
		self.score = score

	def __str__(self):
		"""The name and class of the Profile"""
		out = self.org_name
		if self.label:
			out = out + ": " + self.label
		if self.score:
		  out = out + "\t" + score
		return out

	def __len__(self):
		"""The length of features
		Returns:
		  int
		"""
		return len(self.features)
