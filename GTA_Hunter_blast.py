 # GTA_Hunter pipeline to first blast a query against GTA database to search for GTA homologs before running GTA_Hunter.py
 # Author: Camille Hankel, Taylor Neely
 # Date: 7/5/2016


###############
### IMPORTS ###
###############
import os
import sys 
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO

import argparse
import numpy as np
import time
import pandas as pd 
from bin.filter_extract_fasta_ch import extract_fasta
from bin.Loader import Loader
from bin.Weight import Weight
from bin.Feature import Feature
from bin.SVM import SVM


### CONSTANTS ###
#################
NREPS = 10 # for xval in SVM
PSE_WEIGHT = 0.05 # for pseaac feature
opt_params = {} # contains the parameters for each gene that yielded the best xval results
				# in the format [k-mers,PseAAC,c-value,t-value]
opt_params['2'] = [2,None,1,.01]
opt_params['3'] = [2,None,1,.01]
opt_params['4'] = [2,None,1,.01]
opt_params['5'] = [3,None,10000,.04]
opt_params['6'] = [3,15,1,0]
opt_params['8'] = [2,3,.1,.01]
opt_params['12'] = [4,5,.1,0]
opt_params['13'] = [2,None,1,.01]
opt_params['14'] = [2,None,1,0]
opt_params['15'] = [2,None,1,0]
opt_params["endo"] = [2,10,100,0]
opt_params["rcc171"] = [4,10,.1,0]


def my_main():
	parser = get_args()
	args = parser.parse_args()
	if args.blast:							
		# Initiate blast wrapper
		if not args.queries and not args.folder:
			print ("please specify a query file or a folder of .faa files to blast and classify")
		else:
			if not args.outdir:
				print ("please specify a directory in which to store output files")
			else:
				#Don't use user-specified weighting or training files
				if args.weight:
					print ("Specified weighting files will be ignored, must first determine which GTA gene homologs you have.")
				if args.gta or args.virus:
					print ("Specified training files will be ignored, must first determine which GTA gene homologs you have.")


				if args.folder:
					org = args.folder[0]
					in_folder = org
					out_folder = args.outdir[0]
					for child in os.listdir(in_folder):
						if child[-3:] == "faa":
							name = "_"+child[:-4]
							query = os.path.join(in_folder,child)
							args.queries = [query]
							args.outdir = [out_folder]
							run_wrapper(args,name)

				else:
					run_wrapper(args)

	else: 
		#run GTA_Hunter without doing a BLAST search 
		if not args.gta or not args.virus:
			print ("You must specify the training files if you are running GTA_Hunter without blast")
		else:
			new_hunter(args)

def get_args():
	# Init parser
	parser = argparse.ArgumentParser(description="Gene Classification Using SVM.")

	### Define Args ###

	# Main
	parser.add_argument("-b", "--blast", action="store_true",
		dest="blast", required=False,
		help="Run blast to identify gta homologs.")
	parser.add_argument("-g", "--GTA", type=str, nargs=1,
		dest="gta", required=False,
		help="The .faa or .fna training file for GTA genes.")
	parser.add_argument("-v", "--virus", type=str, nargs=1,
		dest="virus", required=False,
		help="The .faa or .fna training file for viral genes.")
	parser.add_argument("-q", "--queries", type=str, nargs=1,
		dest="queries", required=False,
		help="The .faa or .fna query file to be classified.")
	parser.add_argument("-o", "--outdir", type=str, nargs=1,
		dest="outdir", required=False,
		help="The folder path in which to store output.")
	parser.add_argument("-k", "--kmer", type=int, nargs="?",
		dest="kmer", required=False, const=4, default=None,
		help="The kmer size needed for feature generation (default=4).")
	parser.add_argument("-p", "--pseaac", nargs="?", type=int,
		dest="pseaac", required=False, const=3, default=None,
		help="Expand feature set to include pseudo amino acid composition. Specify lamba (default=3). Weight = 0.05.")
	parser.add_argument("-y", "--physico", action="store_true",
		dest="physico", required=False,
		help="Expand feature set to include physicochemical composition.")
	parser.add_argument("-m", "--min", action="store_true",
		dest="mini", required=False,
		help="Print bare minimum results.")
	parser.add_argument("-O", "--optimal",action="store_true",
		dest="opt", required=False,
		help="Pick the optimal parameters for the GTA gene homolog classification")
	parser.add_argument("-f", "--folder",type=str, nargs=1,
		dest="folder", required=False,
		help="Feed GTA_Hunter a folder from which to run the program on all ")

	# Weight
	parser.add_argument("-W", action="store_true",
		dest="wt", required=False,
		help="Weight training set if desired. Distance files will be supplied automatically.")
	parser.add_argument("-w", "--weight", type=str, nargs=2,
		dest="weight", required=False,
		help="Allows for weighting of training set. Will need to specify the two pairwise distance files needed for weighting (GTA first, then virus).")
	parser.add_argument("-t", "--cluster_type", type=str, nargs=1,
		dest="cluster_type", required=False, default=['farthest'],
		help="Specify 'farthest' or 'nearest' neighbors clustering (default='farthest').")
	parser.add_argument("-d", "--dist", type=float, nargs=1,
		dest="dist", required=False, default=[0.01],
		help="Specify the cutoff distance for clustering in the weighting scheme (default=0.01).")

	# SVM
	parser.add_argument("-c", "--soft_margin", type=float, nargs=1,
		dest="c", required=False, default=[1.0],
		help="The soft margin for the SVM (default=1.0).")
	parser.add_argument("-x", "--xval", type=int, nargs="?",
		dest="xval", required=False, const=5, default=None,
		help="Performs cross validation of training set. Specify folds over 10 repetitions (default=5).")
	parser.add_argument("-e", "--kernel", nargs=2,
		dest="kernel", required=False, default=["linear",0],
		help="Specify kernel to be used and sigma if applicable (i.e. gaussian) (default='linear', 0).")
	parser.add_argument("-s", "--svs", action="store_true",
		dest="svs", required=False,
		help="Show support vectors.")


	return parser

def get_dict():
	#function to create a dictionary with info associated with each gta gene, 
	#so that when there is a blast hit we can identify which gta gene it is a homolog to
	folder = "data/training/gta/"
	folder2 = "data/training/viral/"
	master_dict = {}
	for child in os.listdir(folder):
		if child[-3:] == "faa":
			file = os.path.join(folder,child)
			gene = child.split('_')[0]
			glist = []
			f = open(file)
			content = f.readlines()
			for line in content:
				if line[0] == '>':
					words = line.split()
					ID = words[0][1:]
					glist.append(int(ID))

			master_dict[gene] = glist

	for child in os.listdir(folder2):
		if child[-3:] == "faa":
			file = os.path.join(folder2,child)
			gene = child.split('_')[0]
			f = open(file)
			content = f.readlines()
			for line in content:
				if line[0] == '>':
					words = line.split()
					ID = words[0][1:]
					master_dict[gene].append(int(ID))

	return master_dict

def run_wrapper(args,name=''):
	my_dict = get_dict()

	query_file = args.queries[0]
	out_dir = args.outdir[0]
	#output file (btab) from the blast against GTA db
	blast_out = out_dir +"/blast"+name + ".out"
	#special blast outformat parameters
	outformat = "'6 qseqid sstart send sframe sseqid pident qlen slen length mismatch gapopen qstart qend evalue bitscore'"
	#run blast search using database of viral and GTA training set
	blastp_cline = NcbiblastpCommandline(query=query_file, db="data/GTA_db/GTA_viral.db", evalue=0.001, outfmt=outformat, out=blast_out, num_threads=2,dbsize=10000000)
	blastp_cline
	stdout,stderr = blastp_cline()
	result = open(blast_out)
	lines = result.readlines()
	result.close()
	handle_in = open(query_file)
	#keep track of which sequences did not have a GTA homolog
	no_hit_list = []
	for record in SeqIO.parse(handle_in, "fasta"):
		no_hit_list.append(record.id)
	handle_in.close()

	#continue if blast search had results
	if len(lines) > 0:
		#file to store the results of running GTA_Hunter
		results_file = out_dir + "/results" +name+".out"
		#run script to extract fasta files from blast search
		genes_found = extract_fasta(blast_out,query_file,out_dir,my_dict,name)

		if len(genes_found) > 0:
			results_handle = open(results_file, 'w')
			# for each GTA gene that had homologs, run GTA Hunter
			for gene in genes_found:
				out_faa = out_dir + "/gta_homolog_" +gene+name+".faa"
	
				#get training files corresponding to the gene identified
				gta_file = "data/training/gta/" + gene + "_gta.faa"
				virus_file = "data/training/viral/" + gene + "_viral.faa"
				gta_weight = "data/training/gta/" + gene + "_gta.dist"
				virus_weight = "data/training/viral/" + gene + "_viral.dist"
				homolog_file = out_faa
				args.queries=[homolog_file]
				args.gta = [gta_file]
				args.virus = [virus_file]
				if args.wt:
					args.weight = [gta_weight,virus_weight]
				if args.opt:
					#get the optimal parameters 
					args.kmer = opt_params[gene][0]
					args.pseaac = opt_params[gene][1]
					c_value = opt_params[gene][2]
					args.c = [c_value]
					d_value = opt_params[gene][3]
					args.dist = [d_value]


					print("running with params",args.kmer,"Pseaac",args.pseaac, 
						c_value,d_value)
				#run GTA_Hunter with all the correct files
				print("running GTA Hunter for gene",gene)
				new_hunter(args,results_handle)

			results_handle.close()
			#close original query file and GTA_Hunter results file
			results_handle = open(results_file)
			lines = results_handle.readlines()

			#update the list of sequences that did not have a GTA homolog
			for line in lines:
				words = line.split()
				if words[0] != "Gene":
					if words[0][1:] in no_hit_list:
						no_hit_list.remove(words[0][1:])
			results_handle.close()

	else:
		#if the blast.out file is empty
		print("Sorry, no GTA homologs were found")

	#write out a list of gene IDs of genes that were not homologs to any GTA gene 
	no_hit_file = out_dir+ "/no_homologs" +name+".txt"
	no_hit_handle = open(no_hit_file,'w')
	for nohit in no_hit_list:
		no_hit_handle.write(nohit + '\n')
	no_hit_handle.close()

def new_hunter(args,results_handle=''):

	start = time.time()
	# Print detail
	mini = args.mini
	### Load training set and make features ###
	gta_file = args.gta[0]
	virus_file = args.virus[0]
	# Load profiles
	gta_profs = Loader.load(gta_file, "GTA")
	viral_profs = Loader.load(virus_file, "virus")
	# Make features
	feats = Feature(gta_profs.profiles + viral_profs.profiles)
	if args.kmer:
		kmer_size = args.kmer
		feats.make_kmer_dict(kmer_size)
		feats.kmer_feat()
	if args.pseaac!=None:
		feats.pseaac(lam=int(args.pseaac), weight=PSE_WEIGHT)
	if args.physico:
		feats.physicochem()

	if not args.kmer and args.pseaac==None and not args.physico:
		print("You must specify at least one feature type (-k, -p, -y).")

	else:
		# Weight if needed
		if args.wt:
			# Get distance threshold
			d = args.dist[0]
			# Get cluster type
			cluster_type = args.cluster_type[0]
			# Weight GTA
			pairwiseGTA = Weight.load(args.weight[0])
			GTA_weight = Weight(gta_profs, pairwiseGTA)
			GTA_clusters = GTA_weight.cluster(cluster_type, d)
			GTA_weight.weight(GTA_clusters)
			# Weight Virus
			pairwiseViral = Weight.load(args.weight[1])
			virus_weight = Weight(viral_profs, pairwiseViral)
			virus_clusters = virus_weight.cluster(cluster_type, d)
			virus_weight.weight(virus_clusters)

		# Create SVM
		c = args.c[0]
		kernel = args.kernel[0]
		kernel_var = float(args.kernel[1])

		svm = SVM(gta_profs, viral_profs, c, kernel, kernel_var)

		# Print support vectors
		if args.svs:
			svm.show_svs()

		# Xval	
		if args.xval:
			nfolds = args.xval
			if args.wt:
				result = svm.xval(nfolds, NREPS, pairwiseGTA, pairwiseViral, cluster_type, d)
			else:
				result = svm.xval(nfolds, NREPS)
			if mini:
				print("GTA Correct\tViral Correct")
				print("%.2f\t%.2f" % (result[0], result[1]))
			else:
				print("We correctly classified (on average) %.2f/%d GTA and %.2f/%d Viral genes." 
				% (result[0], len(gta_profs), result[1], len(viral_profs)))
		
		else:
			if args.queries == None:
				print("The query file was not specified. Please declare queries using -q.")
			else: # All good
			# Load test set
				test_profs = Loader.load(args.queries[0])
				# Make features
				if args.kmer:
					feats.kmer_feat(test_profs)
				if args.pseaac:
					feats.pseaac(lam=int(args.pseaac), weight=PSE_WEIGHT, profiles=test_profs)
				if args.physico:
					feats.physicochem(profiles=test_profs)
				# Classify
				svm.predict(test_profs)
				# The gta gene being classified
				gene = args.gta[0].split('/')[3].split('_')[0]
				# Print results
				if mini:
					print("Gene\t\tClass")
					if args.blast:
						results_handle.write("Gene\t\tClass\n")
					for profile in test_profs:
						print(">%s\t%s" % (profile.name, profile.label))
						if args.blast:
							results_handle.write(">%s\t%s\n" % (profile.name, profile.label))
				else:
					print("%-*s%-*s%-*s%12s" % (55, "Gene", 10, "Score", 5, "Classification","GTA Gene"))
					if args.blast:
						results_handle.write("%-*s%-*s%-*s%12s\n" % (55, "Gene", 10, "Score", 5, "Classification", "GTA Gene"))
					for profile in test_profs:
						print(">%-*s%-*f%-*s%14s" % (55, profile.org_name, 10, profile.score, 5, profile.label, gene))
						if args.blast:
							results_handle.write(">%-*s%-*f%-*s%14s\n" % (55, profile.org_name, 10, profile.score, 5, profile.label, gene))

	end = time.time()
	total = (end - start)
	print ("time to run:", total)

my_main()