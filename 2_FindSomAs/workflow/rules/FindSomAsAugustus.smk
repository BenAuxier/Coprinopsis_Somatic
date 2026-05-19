### FindSomAs: Find homologs of the somA locus in Coprinopsis cinerea and friends
#############################################################################
#############################################################################
# ==================================================
# S. Lorena Ament-Velasquez
# 2024-12-10
# ==================================================

# -------------------------------------------------
# Variables from configuration file
GenomesIDs = config["GenomesIDs"]

# path2genomes = config["path2genomes"]
path2hmmers = config["path2hmmers"]

# Variables
FLANKs = config["FLANKs"]

# Scripts
hmmer2bed = config["hmmer2bed"]
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: 
# ----------

## ---- Prepare files -----
# rule prepare_hmmeroutput:
# 	output:
# 		"data/HMMERresults/{prot}_hits_{sample}.txt"
# 	shell:
# 		"cp {path2hmmers}/{wildcards.prot}_hits/*_hits*{wildcards.sample}* {output}"

## -------------------------
rule hmmer2bed:
	input:
		gff = "data/Annotations/{sample}.gff",
		hmmer = path2hmmers + "/{prot}_hits_{sample}.txt",
	output:
		"workflow/intermediate/HMMERhits/{prot}_hits_{sample}.bed",
	params:
		flanks = FLANKs
	shell:
		"python {hmmer2bed} {input.gff} {input.hmmer} -f {params.flanks} -a > {output}" # Augustus specific

rule BEDtools:
	input:
		nacht = "workflow/intermediate/HMMERhits/nacht_hits_{sample}.bed",
		kinase = "workflow/intermediate/HMMERhits/kinase_hits_{sample}.bed",
	output:
		"workflow/intermediate/beds/{sample}.overlaps"
	shell:
		"bedtools intersect -a {input.nacht} -b {input.kinase} -wo > {output}"

def remove_overlap(ranges):
	""" Simplify a list of ranges; I got it from https://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap """
	result = []
	current_start = -1
	current_stop = -1 

	for start, stop in sorted(ranges):
		if start > current_stop:
			# this segment starts after the last segment stops
			# just add a new segment
			result.append( (start, stop) )
			current_start, current_stop = start, stop
		else:
			# current_start already guaranteed to be lower
			current_stop = max(current_stop, stop)
			# segments overlap, replace
			result[-1] = (current_start, current_stop) # SLAV: I modified this to update the stop too. (otherwise small ranges contained in a previous larger range will make the current stop smaller than it should be)
	return(result)

rule simplify_bed:
	input:
		bed = "workflow/intermediate/beds/{sample}.overlaps",
	output:
		bed = "results/beds/{sample}_somA_candidates.bed"
	run:
		tabopen = open(input.bed, 'r')
		tabs = [line.rstrip("\n").split("\t") for line in tabopen]

		chunks = {} # Collect the ranges of all regions that have nachts and kinases nearby each other
		for tab in tabs:
			ctg = tab[0]
			start = int(tab[1])
			end = int(tab[2])

			if ctg in chunks.keys():
				chunks[ctg].append((start, end))
			else:
				chunks[ctg] = [(start, end)]

		ofile = open(output.bed, 'w')

		# Reduce overlapping ranges and print them
		for ctg in chunks.keys():
			reducedranges = remove_overlap(chunks[ctg])
			for rg in reducedranges:
				ofile.write(f"{ctg}\t{rg[0]}\t{rg[1]}\n")


rule get_GFFSlicer:
	output:
		"workflow/scripts/GFFSlicer.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/refs/heads/master/GenomeAnnotation/GFFSlicer.py"


rule GFFSlicer:
	""" Subset the gff file based on the candidate regions """
	input:
		GFFSlicer = "workflow/scripts/GFFSlicer.py",
		bed = "results/beds/{sample}_somA_candidates.bed",
		gff = "data/Annotations/{sample}.gff",
	output:
		gff = temp("results/gffs/{sample}_somA_candidates.gff-raw"),
		report = temp("workflow/temp/{sample}_somA_candidates-raw.txt")
	run:
		tabs = [line.rstrip("\n").split("\t") for line in open(input.bed, 'r')]

		reportfile = open(output.report, 'w')
		reportfile.write(f"{wildcards.sample}\t{len(tabs)}\n")

		ofile = open(output.gff, 'w')
		ofile.write("##gff-version 3\n")
		ofile.close()

		for tab in tabs:
			ctg, start, end = tab

			cmd = f"python {input.GFFSlicer} {input.gff} {start} {end} --contig {ctg} --outputname {wildcards.sample}_{ctg}.gff --keepcoords --nocomments"
			shell(cmd)
			shell(f"cat {wildcards.sample}_{ctg}.gff >> {output.gff}")
			shell("rm {wildcards.sample}_{ctg}.gff")

rule pretty_gff:
	""" Add notes and colors to the genes so they are easier to see in IGV """
	input:
		gff = "results/gffs/{sample}_somA_candidates.gff-raw",
		nacht = "workflow/intermediate/HMMERhits/nacht_hits_{sample}.bed",
		kinase = "workflow/intermediate/HMMERhits/kinase_hits_{sample}.bed",
		report = "workflow/temp/{sample}_somA_candidates-raw.txt"
	output:
		gff = "results/gffs/{sample}_somA_candidates.gff",
		report = "workflow/temp/{sample}_somA_candidates.txt"
	run:
		# Read the list of genes identified by HMMER
		tabs_nacht = [line.rstrip("\n").split("\t") for line in open(input.nacht, 'r')]
		tabs_kinase = [line.rstrip("\n").split("\t") for line in open(input.kinase, 'r')]

		nachtids = [tab[3] for tab in tabs_nacht] # Use the ID of the transcript, Augustus-specific
		kinaseids = [tab[3] for tab in tabs_kinase]

		# Read the report of how many loci
		loci = [line.rstrip("\n").split("\t") for line in open(input.report, 'r')][0][1]

		list_nachts = []
		list_kinases = []
		with open(input.gff, 'r') as file, open(output.gff, 'w') as ofile: # Start a gff to put all the genes
			for line in file:
				if ('\tgene\t' in line) or ('\ttranscript\t' in line): # Only modify the genes and mRNA features
					nakedline = line.rstrip("\n")
					tabs_gene = nakedline.split("\t")
					attributes = tabs_gene[8].split(";")

					for att in attributes:
						if 'ID=' in att:
							locus_tag = att.replace('ID=', '')
					if '.t1' not in locus_tag: locus_tag += '.t1'

					if locus_tag in nachtids:
						ofile.write(nakedline + ';Alias=NACHT;color="#ffd700";\n')
						list_nachts.append(locus_tag)
					elif locus_tag in kinaseids:
						ofile.write(nakedline + ';Alias=Kinase;color="#ff4500";\n')
						list_kinases.append(locus_tag)
					else:
						ofile.write(line)
				else:
					ofile.write(line)

		# Remove redundancy and sort
		list_nachts = list(set(list_nachts))
		list_kinases = list(set(list_kinases))
		list_nachts.sort()
		list_kinases.sort()

		# Report
		with open(output.report, 'w') as reportfile:
			reportfile.write(f"{wildcards.sample}\t{loci}\t{len(list_nachts)}\t{len(list_kinases)}\t{','.join(list_nachts)}\t{','.join(list_kinases)}\n")

rule report_counts:
	""" Make a report of the results """ 
	input:
		expand("workflow/temp/{sample}_somA_candidates.txt", sample = GenomesIDs)
	output:
		"results/Loci_counts_Augustus.txt"
	shell:
		"echo 'Strain\tN_candidate_loci\tN_NACHTs\tN_Kinases\tNACHTS_locustags\tKinases_locustags' > {output}; "
		"cat {input} >> {output}"

# ---- Get target genes into alignments to make trees ----

rule get_gffutils2fasta:
	output:
		"workflow/scripts/gffutils2fasta.py"
	shell:
		"wget -O {output} https://raw.githubusercontent.com/SLAment/Genomics/refs/heads/master/GenomeAnnotation/gffutils2fasta.py"


rule make_prot_fastas:
	input:
		gff = "data/Annotations/{sample}.gff",
		genome = "data/Genomes/{sample}.fa",
		report = "workflow/temp/{sample}_somA_candidates.txt",
		gffutils2fasta = "workflow/scripts/gffutils2fasta.py"
	output:
		nachts = "workflow/intermediate/alignments/{sample}_somA_candidates_nacht.faa",
		kinases = "workflow/intermediate/alignments/{sample}_somA_candidates_kinase.faa",
	run:
		tabs = [line.rstrip("\n").split("\t") for line in open(input.report, 'r')][0]

		try:
			NACHTs = tabs[4]
			Kinases = tabs[5]
		except:
			NACHTs = []
			Kinases = []

		if len(NACHTs) > 0:
			# Extract the protein sequence of those genes
			cmd = f"python {input.gffutils2fasta} {input.genome} {input.gff} --type CDS --join -p --specificgene {NACHTs} --onlyids --mRNAids --output {output.nachts}"
			shell(cmd)
			shell(f"sed -i 's/>/>{wildcards.sample}_/' {output.nachts}")

			cmd = f"python {input.gffutils2fasta} {input.genome} {input.gff} --type CDS --join -p --specificgene {Kinases} --onlyids --mRNAids --output {output.kinases}"
			shell(cmd)
			shell(f"sed -i 's/>/>{wildcards.sample}_/' {output.kinases}")
		else:
			shell("touch {output.nachts}")
			shell("touch {output.kinases}")

rule cat_prot_fastas:
	input:
		expand("workflow/intermediate/alignments/{sample}_somA_candidates_{{prot}}.faa", sample = GenomesIDs)
	output:
		"results/fastas/somA_candidates_{prot}.faa"
	shell:
		"cat {input} | sed 's/*//g' > {output}; "

