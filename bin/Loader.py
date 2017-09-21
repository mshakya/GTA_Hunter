"""
File name: Loader.py
Date created: 10/22/2015
Date last modified: 05/31/2016
Python version: 3.5.1
Description: A class containing a list
	of profiles, loaded in by file.
	Handle all file checking.
"""

###############
### IMPORTS ###
###############
from Profile import Profile
import os
import argparse
import numpy as np

############
### CODE ###
############
class Loader:
	"""Contains all information pertaining to a given gene."""

	def __init__(self, profiles, profileDict, label):
		""" Do not call directly, use the load 
			method instead.
			Example:
				sequences = Loader.load("file.fa", "GTA")
		"""
		self.profiles = profiles
		self.profileDict = profileDict
		self.label = label

	@staticmethod
	def load(filename, label=None):
		"""Loads given file into list of Profiles
		Args:
		  file (string): the name of the file to be loaded
		  label (string): the assigned label for the file, 
			used in training SVM (optional)
		Returns:
			A Loader Object with the loaded profiles
		Example:
			profiles = Loader.load('genes/gene5.fna', label='GTA')
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
		profiles = []
		profileDict = {}
		# Parse file
		isFirst = True
		for line in open(filename, 'r'):
			line = line.strip()
			# Get info if not the start of an element
			if line[0] != '>':
				if target == "features": # csv
					data = line.split(',')
					# convert to int/float
					for i in range(len(data)):
						data[i] = float(data[i])
				elif target == "weight": # float
					data = float(line)
				else: # multiline string
					data += line
			# Start of an element, initialize vars
			else:
				# update data and save profile of previous element
				if not isFirst:
					setattr(profile, target, data)
					profiles.append(profile)
					profileDict[name] = profile
				# init new profile info
				name = line[1:].split()[0]
				org_name = line[1:]
				# if '[' in line and ']' in line:
				# 	org_name = line[line.index('[')+1:line.index(']')]
				# else:
				# 	org_name = name
				if target == "prot_seq" or target == "dna_seq":
					data = ''
				# create profile
				profile = Profile(name, org_name, label)
				# Clear isFirst
				isFirst = False

		# update data and save profile for last one
		setattr(profile, target, data)
		profiles.append(profile)
		profileDict[name] = profile

		return Loader(profiles, profileDict, label)


	def __getitem__(self, pos):
		"""The Profile at given postion
		Returns:
		  Object: Profile
		"""
		if type(pos) == int:
			return self.profiles[pos]
		elif type(pos) == str:
			return self.profileDict[pos]

	def __len__(self):
		"""The number of profiles loaded
		Returns:
		  int
		"""
		return len(self.profiles)

	def __str__(self):
		"""The class of the profiles"""
		return self.label

	def write(self, profiles, datatype, outfile):
		"""Writes data to target outfile"""
		out = open(outfile, 'w')
		for profile in profiles:
			out.write('>' + profile.org_name + '\n')
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

if __name__ == '__main__':
	# Define arg parser
	parser = argparse.ArgumentParser(description="Test uniquing.")

	### Define Args ###

	# Main
	parser.add_argument("-g", "--GTA", type=str, nargs=1,
		dest="gta", required=True,
		help="The .faa or .fna training file for GTA genes.")
	parser.add_argument("-v", "--virus", type=str, nargs=1,
		dest="virus", required=True,
		help="The .faa or .fna training file for viral genes.")
	parser.add_argument("-q", "--queries", type=str, nargs=1,
		dest="queries", required=True,
		help="The .faa or .fna query file to be purged.")
	parser.add_argument("-o", "--out", type=str, nargs=1,
		dest="outfile", required=True,
		help="The .faa or .fna out file to be written.")

	args = parser.parse_args()
	
	### Load training set and testset ###
	gta_file = args.gta[0]
	virus_file = args.virus[0]
	test_file = args.queries[0]
	out_file = args.outfile[0]
	# Load profiles
	gta_profs = Loader.load(gta_file, "GTA")
	viral_profs = Loader.load(virus_file, "virus")
	test_profs = Loader.load(test_file, "test")
	
	# Add unique to list
	uniques = []
	for profile in test_profs:
		if profile.name in gta_profs.profileDict or profile.name in viral_profs.profileDict:
			# print("I'm not unique: %s" % profile)
			pass
		else:
			uniques.append(profile)
	# Write out uniques from test set
	test_profs.write(uniques, "prot_seq", out_file)