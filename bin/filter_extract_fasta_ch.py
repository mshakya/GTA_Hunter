#!/usr/bin/env python
#File created on April 3 2014

"""the script formats the blastp results
"""

__author__ = "Migun Shakya", "Camille Hankel"
__email__ = "microbeatic@gmail.com"
__version__ = "0.1"
__license__ = "The MIT License (MIT)"


#--- standard library imports
import argparse
import sys
import os

#--- third-party imports
import pandas as pd
import numpy as np
#--- project specific imports


def cmdline_parser():
    """
    creates an argparse instance
    """
    parser = argparse.ArgumentParser(description="""filters the blast results
                                     the program also filters out tblastn hits with less than 0.6 alignment
                                     length, and performs furhter filtering on genes with hits in same region""")
    parser.add_argument("--blast_hits",
                        help="""blast results""",
                        dest="blast_hits",
                        required=True)
    parser.add_argument("--header",
                        help="""header corresponding to results of blast""",
                        dest="header",
                        required=True)
    parser.add_argument("--fasta",
                        help="source fasta",
                        required=True)
    parser.add_argument("--out_fasta",
                        help=""" output fasta file""",
                        dest="out_fasta",
                        required=True)
    parser.add_argument("--db_source", help="""parse fasta from ncbi or img
                        """, choices=['img', 'ncbi', 'phantome'], default='ncbi')
    parser.add_argument("-aln_len", help="""proportion of query
                        alignment used by""",
                        type=float,
                        default=0.6)
    parser.add_argument("-slen", help="""proportion of subject length
                        used in the alignment""",
                        type=float,
                        default=0.6)
    return parser


def extract_fasta(blast_out,src_fasta,out_dir,my_dict,name):
    """
    The main function
    """
    # parser = cmdline_parser()
    # args = parser.parse_args()
    #check if the files are empty or exist first
    try:
        if os.stat(blast_out).st_size > 0:
            print ("Processing...")
        else:
            print ("Empty URL file ... exiting")
            sys.exit()
    except OSError:
        print ("URL file missing ... exiting")
        sys.exit()

    #1 read blast results into pandas dataframe
    blast_hits = pd.read_csv(blast_out,
                             sep='\t',
                             header=None,
                             )

    #2 add header to the dataframe
    header_list = ["qseqid","sstart", "send", "sframe", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend", "evalue", "bitscore"]
    print (header_list)
    blast_hits.columns = header_list

    #3 add columns with proportion of query that aligned with the subject
    # and the proportion of subject that aligned with the query
    blast_hits['proportion'] = blast_hits.apply(
        lambda row: ((float(row['qend'])-float(row['qstart']))/float(row['qlen'])), axis=1)
    blast_hits['sprop'] = blast_hits.apply(
        lambda row: ((float(row['send'])-float(row['sstart']))/float(row['slen'])), axis=1)

    print (blast_hits.head(25))

    #4 filter data frame based on the proportion
    #remove rows with less than arg.aln_len of query used in the alignment
    blast_hits_v1 = blast_hits[(blast_hits['proportion'] > .4)]
    blast_hits_v2 = blast_hits_v1[(blast_hits['sprop'] > .4)]
 
    print (blast_hits_v1.head(25))
    print (blast_hits_v2.head(25))

    # blast_hits_v3 = blast_hits_v2.drop_duplicates('qseqid')
    blast_hits_v3 = blast_hits_v2
    print (blast_hits_v3.head(25))

    found_genes = []
    if blast_hits_v3.empty:
        print("sorry no GTA homologs")
    else:
        for key in my_dict.keys():
            blast_hits_gene = blast_hits_v3.loc[blast_hits_v3['sseqid'].isin(my_dict[key])]
            if not blast_hits_gene.empty:
                gi_list = blast_hits_gene['qseqid'].tolist()
                gi_list = [str(i) for i in gi_list]
                out_fasta = out_dir + "/gta_homolog_" +key+name+".faa"
                found_genes.append(key)
                phantome2fasta(src_fasta, out_fasta, gi_list)

    print (found_genes)
    return (found_genes)

def img2fasta(source_fasta, result_fasta, gi_list):
    """ extracts sequences from fasta A with
    gi numbers in the given list B
    """
    from Bio import SeqIO
    out_fasta = open(result_fasta, 'w')
    open_fasta = SeqIO.parse(open(source_fasta), 'fasta')
    for seq in open_fasta:
        if seq.id in gi_list:
            if seq.seq[0] == 'M':
                print (seq.id)
                SeqIO.write(seq, out_fasta, 'fasta')
    out_fasta.close()


def gi2fasta(source_fasta, result_fasta, gi_list):
    """ extracts sequences from fasta A with
    gi numbers in the given list B
    """
    from Bio import SeqIO
    out_fasta = open(result_fasta, 'w')
    open_fasta = SeqIO.parse(open(source_fasta), 'fasta')
    for seq in open_fasta:
        if seq.id.split('|')[1] in gi_list:
            SeqIO.write(seq, out_fasta, 'fasta')
    out_fasta.close()


def phantome2fasta(source_fasta, result_fasta, gi_list):
    """ extracts sequences from fasta A with
    gi numbers in the given list B
    """
    from Bio import SeqIO
    out_fasta = open(result_fasta, 'w')
    open_fasta = SeqIO.parse(open(source_fasta), 'fasta')
    for seq in open_fasta:
        if seq.id in gi_list:
            SeqIO.write(seq, out_fasta, 'fasta')
    out_fasta.close()


def get_dict():
    #function to create a dictionary with info associated with each gta gene, 
    #so that when there is a blast hit we can identify which gta gene it is a homolog to
    folder = "data/training/gta/"
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

    return master_dict


if __name__ == '__main__':

    my_dict = get_dict()
    extract_fasta("heatmaps/blast3.out","test45.faa","heatmaps",my_dict)

    if sys.version_info < (2, 7) or sys.version_info > (2, 8):
        LOG.info("only tested Python 2.7 so far")

