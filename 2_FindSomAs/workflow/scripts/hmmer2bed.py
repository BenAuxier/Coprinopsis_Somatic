#!/usr/bin/env python
# encoding: utf-8

# ================== hmmer2bed =================
# Take a hmmer output and the GenBank annotation of an assembly genome, and
# output a BED file with the gene coordinates of the proteins found by hmmer

# ==================================================
# Sandra Lorena Ament Velasquez
# 2024/12/10
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
import sys  # To exit the script
import gffutils
import re
# ------------------------------------------------------

version = 1.2
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Make a BED file with the gene coordiantes of the proteins found by hmmer *", epilog="")  # Create the object using class argparse

# Add options
parser.add_argument('GFF', help="GFF3 file from GenBank")
parser.add_argument('HMMER', help="Output file of hmmer")

# Options
parser.add_argument("--evalue", "-e", help="Minimum e-value of hmmer hit to be considered (default 0.00001)", type=float, default=0.00001)
parser.add_argument("--flanks", "-f", help="Add a given number of bp to the right and left flanks of each gene (default 0)", type=int, default=0)
parser.add_argument("--augustus", "-a", help="The GFF file comes from Augustus", default=False, action='store_true')

# extras
# parser.add_argument('--output', '-o', help="Name of output file (default is to append 'ed' to the input name)", default=None)
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)


try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	GFFopen = open(args.GFF, 'r')
	HMMERopen = open(args.HMMER, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------


# ---------------------------------
# Make database
# ---------------------------------
# This will parse the file, infer the relationships among the features in the file, and store the features and relationships
# See https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html

dbfnchoice = ':memory:'

# http://daler.github.io/gffutils/database-ids.html
id_spec={"gene": ["ID", "Name"], 
	"mRNA": ["ID", "transcript_id"], 
	"transcript": ["ID", "transcript_id"], 
	"rRNA": ["ID", "Name"],
	"tRNA": ["ID", "Name"]} 

db = gffutils.create_db(data = args.GFF, 
	keep_order = True,
	dbfn = dbfnchoice,
	force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	verbose = False,
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 

# ---------------------------------

# --- Read the hmmer output ---
hmmerprots = {} # make a dictionary to hold the protein_ids and the e-value of the full protein according to HMMER
with HMMERopen as file:
	for line in file:
		if not line.startswith("#"):  # Skip lines that start with the '#' symbol
			line_with_tabs = re.sub(r"\s+", "\t", line) # replace one or more white spaces with a tab
			tabs = line_with_tabs.rstrip("\n").split("\t")

			protein_id = tabs[0]
			proteval = float(tabs[4]) #e-value of the full sequence

			if proteval < args.evalue:
				hmmerprots[protein_id] = proteval

# --- Read the GFF  ---
# sys.stdout.write("contig\tstart\tend\tprotein_id\tproduct\tevalue\n") # head
# GeneProtDic = {}
for gene in db.features_of_type('gene'):
	geneID = gene['ID'][0]
	contig = gene.chrom
	start = int(gene.start) - args.flanks
	
	if start < 0: start = 0 # sometimes a gene can be at the start of the contig and create negative numbers if the flanks are modified
	end = gene.end + args.flanks

	if args.augustus: # Assume it's an Augustus gff
		protname = geneID + '.t1'

		if protname in hmmerprots.keys():
			sys.stdout.write(f"{contig}\t{start}\t{end}\t{protname}\t.\t{hmmerprots[protname]}\t.\n")

	else: # Assume it's a GenBank gff
		for child in db.children(gene, featuretype='CDS', order_by='start'):
			protname = child['Name'][0]
			product = child['product'][0]
			locus_tag = child['locus_tag'][0]
			break

		# sys.stdout.write(line)
		# GeneProtDic[geneID] = (protname, contig, start, end, product)
		if protname in hmmerprots.keys():
			# print(contig, start, end, protname) # simple BED file
			sys.stdout.write(f"{contig}\t{start}\t{end}\t{protname}\t{product}\t{hmmerprots[protname]}\t{locus_tag}\n")

