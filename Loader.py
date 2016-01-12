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

	def __init__(self, profiles, label):
		""" Do not call directly, use the load 
			method instead.
			Example:
				sequences = Loader.load("file.fa", "GTA")
		"""
		self.profiles = profiles
		self.label = label

	@staticmethod
	def load(filename, label=None):
		"""Loads given file into list of Profiles
		Args:
		  file (string): the name of the file to be loaded
		  label (string): the assigned label for the file, 
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
		isFirst = True
		for line in open(filename, 'r'):
			line = line.strip()
			if line[0] == '>':
				name = line[1:]
				if target == "prot_seq" or target == "dna_seq":
					data = ''
				if not isFirst: # update data and save profile
					setattr(profile, target, data)
					profile[name] = profile
				# create profile
				profile = Profile(name, label)
				isFirst = False
			else:
				if target == "features": # csv
					data = line.split(',')
					# convert to int/float
					for i in range(len(data)):
						data[i] = int(data[i])
						# except ValueError:
						# 	data[i] = float(data[i])
				elif target == "weight": # float
					data = float(line)
				else: # multiline string
					data += line

		# update data and save profile for last one
		setattr(profile, target, data)
		profile[name] = profile

		return Loader(profiles, label)


	def __getitem__(self, pos):
		"""The Profile at given postion
		Returns:
		  Object: Profile
		"""
		return self.profiles[pos]

	def __len__(self):
		"""The number of profiles loaded
		Returns:
		  int
		"""
		return len(self.profiles)

	def __str__(self):
		"""The class of the profiles"""
		return self.label

	def write(self, datatype, outfile):
		"""Writes data to target outfile"""
		out = open(outfile, 'w')
		for profile in profiles:
			out.write('>' + profile.name + '\n')
			if datatype == "prot_seq" or datatype == "dna_seq":
				data = getattr(profile, datatype)
				for i in range(0, len(data), 70):
					out.write(data[i:i+70] + '\n')
			elif datatype == "features":
				data = getattr(profile, datatype)
				for i in range(len(data)-1):
					out.write(data[i] + ',')
				out.write(data[len(data)-1] + '\n')
			elif datatype == "weight":
				data = getattr(profile, datatype)
				out.write(data + '\n')
			else:
				raise Exception("The datatype %r is not supported."%(datatype))
				out.close()
				return

		out.close()
		return

	

