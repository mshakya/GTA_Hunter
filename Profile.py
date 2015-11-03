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

    def __init__(self, name, classification=None, dna_seq=None, prot_seq=None, weight=None, features=None):
        """A Profile with at least one initialized variable
        Args:
          name (string): name of the gene
          classification (string): classification of gene (optional)
          dna_seq (string): DNA sequence from FASTA (optional)
          prot_seq (string): AA seqeuence from FASTA (optional)
          weight (float): given weight when training with profile (optional)
          features: a sequence of features derived from Feature (optional)
        """
        self.name = name
        self.classification = classification
        self.dna_seq = dna_seq
        self.prot_seq = prot_seq
        self.weight = weight
        self.features = features

    def __str__(self):
        """The name and class of the Profile"""
        out = self.name
        if self.classification:
            out = out + ": " + self.classification
        return out
