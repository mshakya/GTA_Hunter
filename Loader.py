"""
File name: Loader.py
Date created: 22/10/2015
Date last modified: 22/10/2015
Python version: 2.7
Description: A class containing a list
	of profiles, loaded in by file.
	Handle all file checking.
"""

###############
### IMPORTS ###
###############
from Profile import Profile
import os
import numpy as np

############
### CODE ###
############
class Loader:
    """Contains all information pertaining to a given gene."""

    def __init__(self, profiles, classification):
        """ Do not call directly, use the load 
            method instead.
            Example:
                sequences = Loader.load("file.fa", "GTA")
        """
        self.profiles = profiles
        self.classification = classification

    @staticmethod
    def load(filename, classification=None):
        """Loads given file into list of Profiles
        Args:
          file (string): the name of the file to be loaded
          classification (string): the assigned classification for the file, 
            used in training SVM (optional)
        Returns:
            A Loader Object with the 
        Example:
        """

        # Check extensions (.fna, .faa, .w, fe, and .km accepted)
        if filename.endswith(".fna"):
            target = "dna_seq"
        elif filename.endswith(".faa"):
            target = "prot_seq"
        elif filename.endswith(".w"):
            target = "weight"
        elif filename.endswith(".fe"):
            target = "features"
        elif filename.endswith(".km"):
            raise ValueError("The load function does not directly handle .km files. Use expand function instead.")
        else:
            raise ValueError("The file type was not recognized, please use of the the following formats: .fna, .faa, .w, .fe.")

        # Check if file exists
        if not os.path.isfile(filename):
            raise Exception("The file %r does not exist."%(filename))
        # Check if file can be read
        if not os.access(filename, os.F_OK):
            raise Exception("The file %r cannot be read."%(filename))

        # Create profile dictionary
        profiles = {}
        # Parse file
        for line in open(filename, 'r'):
            print line
            # getattr(profile, target)
        #     if line[0:4] == 'ATOM' and line[13:15] == 'CA':
        #         resi = int(line[23:26])
        #         if resi in resis:
        #             print 'duplicate for resi',resi,'-- ignoring'
        #         else:
        #             resis.add(resi)
        #             cas.append(np.array([float(line[31:39]), float(line[39:47]), float(line[47:54])]))
        # return Loader(profiles, classification)


    def __getitem__(self, pos):
        """The Profile at given postion
        Returns:
          Object: Profile
        """
        return self.profiles[pos]

    def __len__(self):
        """The width of the motif
        Returns:
          int
        """
        return len(self.profiles)

    def __str__(self):
        """The name of the Profile"""
        return self.name