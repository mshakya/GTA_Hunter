"""
File name: GTA_Tool.py
Date created: 22/10/2015
Date last modified: 4/12/2015
Python version: 2.7
Description: Combines Loader.py,
	Profile.py, SVM.py, Features.py,
	Weight.py, and Filter.py to be able
	to characterize and classify two different
	types of genes (GTA and phage).
"""

###############
### IMPORTS ###
###############
from Loader import Loader
from Filter import Filter
from Weight import Weight
from Feature import Feature
from SVM import SVM

############
### META ###
############
__author__ = "Taylor Neely"
__version__ = "0.0.1"
__email__ = "tneely@dartmouth.edu"
__status__ = "Development"

############
### CODE ###
############