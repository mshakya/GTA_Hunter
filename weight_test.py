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
from Weight import Weight
from Feature import Feature
from SVM import SVM
import matplotlib.pyplot as plt

############
### CODE ###
############
if __name__ == '__main__':
	
	# Weight Test
	c1 = Loader.load("test/class1.fe")
	c2 = Loader.load("test/class2.fe")
	c1reps = Loader.load("test/class1_reps.fe")
	print(c1reps[26])
	# Set weights
	for p1, p2 in zip(c1,c2):
		p1.weight = 1
		p2.weight = 1
	
	# Make SVM
	svm = SVM(c1, c2, 10.0, "linear", sigma=None, xval=False)
	print("Alphas:",svm.predictor.alphas,
		"\nSupport vectors (%d):\n" % len(svm.predictor.svs),
		svm.predictor.svs,"\nSupport vector labels:",svm.predictor.svsLabels)
	svsx = [sv[0] for sv in svm.predictor.svs]
	svsy = [sv[1] for sv in svm.predictor.svs]

	# Visualize
	c1x = [profile.features[0] for profile in c1]
	c1y = [profile.features[1] for profile in c1]
	c2x = [profile.features[0] for profile in c2]
	c2y = [profile.features[1] for profile in c2]
	fig = plt.figure(4)
	ax1 = fig.add_subplot(111)
	ax1.scatter(c1x, c1y, c='b', marker='s', label='class1')
	ax1.scatter(c2x, c2y, c='r', marker='s', label='class2')
	plt.legend(loc='upper left')
	plt.title('Original plot')

	# Visualize
	c1x = [profile.features[0] for profile in c1]
	c1y = [profile.features[1] for profile in c1]
	c2x = [profile.features[0] for profile in c2]
	c2y = [profile.features[1] for profile in c2]
	fig = plt.figure(1)
	ax1 = fig.add_subplot(111)
	ax1.scatter(c1x, c1y, c='b', marker='s', label='class1')
	ax1.scatter(c2x, c2y, c='r', marker='s', label='class2')
	ax1.scatter(svsx, svsy, c='y', marker='o', label='Svs')
	plt.legend(loc='upper left')
	plt.title('One c1 outlier, no weights (1)')
	
	for p in c1reps:
		p.weight = 1

	# Make SVM
	svm = SVM(c1reps, c2, 10.0, "linear", sigma=None, xval=False)
	print("Alphas:",svm.predictor.alphas,
		"\nSupport vectors (%d):\n" % len(svm.predictor.svs),
		svm.predictor.svs,"\nSupport vector labels:",svm.predictor.svsLabels)
	svsx = [sv[0] for sv in svm.predictor.svs]
	svsy = [sv[1] for sv in svm.predictor.svs]

	# Visualize
	c1x = [profile.features[0] for profile in c1reps]
	c1y = [profile.features[1] for profile in c1reps]
	c2x = [profile.features[0] for profile in c2]
	c2y = [profile.features[1] for profile in c2]
	fig = plt.figure(2)
	ax1 = fig.add_subplot(111)
	ax1.scatter(c1x, c1y, c='b', marker='s', label='class1')
	ax1.scatter(c2x, c2y, c='r', marker='s', label='class2')
	ax1.scatter(svsx, svsy, c='y', marker='o', label='Svs')
	plt.legend(loc='upper left')
	plt.title('6 c1 outliers, no weights (1)')

	# Set weights
	for p in range(25,31):
		c1reps[p].weight = 1/6

	# Make SVM
	svm = SVM(c1reps, c2, 10.0, "linear", sigma=None, xval=False)
	print("Alphas:",svm.predictor.alphas,
		"\nSupport vectors (%d):\n" % len(svm.predictor.svs),
		svm.predictor.svs,"\nSupport vector labels:",svm.predictor.svsLabels)
	svsx = [sv[0] for sv in svm.predictor.svs]
	svsy = [sv[1] for sv in svm.predictor.svs]

	# Visualize
	c1x = [profile.features[0] for profile in c1reps]
	c1y = [profile.features[1] for profile in c1reps]
	c2x = [profile.features[0] for profile in c2]
	c2y = [profile.features[1] for profile in c2]
	fig = plt.figure(3)
	ax1 = fig.add_subplot(111)
	ax1.scatter(c1x, c1y, c='b', marker='s', label='class1')
	ax1.scatter(c2x, c2y, c='r', marker='s', label='class2')
	ax1.scatter(svsx, svsy, c='y', marker='o', label='Svs')
	plt.legend(loc='upper left')
	plt.title('6 c1 outliers, weighted (1/6)')
	plt.show()