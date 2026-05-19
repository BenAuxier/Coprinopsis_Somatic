# Plotting counts of candidate somA loci across the phylogeny

## Making a phylogeny of *Coprinopsis cinerea* and related species

First I got the sequences for all species for 3 genes, using the homolog of *C. cinerea*: the classic phylogenetics marker *rpb2* (Genbank accession number XM_001837441.2), as well as *mcm7* (XM_001833124.2) and *tsr1* (XP_001830678.2), two genes with good phylogenetic signal for fungi ([Aguileta et al., 2008](https://academic.oup.com/sysbio/article-abstract/57/4/613/1632046)). 

I performed tBLASTn searches of these *C. cinerea* references on each species genome with the script [`query2haplotype.py`](https://github.com/SLAment/Genomics/blob/master/BLAST/query2haplotype.py). I ran it like so:

	$ python query2haplotype.py {input.genome} {input.refgene} --task tblastn --haplo --minhaplo 200 --identity 60 --extrabp 1600 | sed "s/>/>{sampleID} /" >> {output.fas}

I aligned each gene with MAFFT v7.526 ([Katoh and Standley, 2013](https://doi.org/10.1093/molbev/mst010)) followed by manual curation using the *C. cinerea* CDS sequences as guides to deal with the very divergent introns. 

The gene alignments were concatenated and given to IQ-TREE v2.3.6 ([Minh et al., 2020](https://doi.org/10.1093/molbev/msaa015)):

	$ iqtree2 -s RPB2+MCM7+TSR1.fa -m MFP -seed 1234 -b 100 -nt 6 -bnni -pre RPB2+MCM7+TSR1

This produced a maximum likelihood phylogeny in newick format.

The alignment in nexus and fasta format, as well as the tree are found in the `data` directory. The nexus format has definition of where the 3 genes start and end. You can visualize it with [SeaView](https://doua.prabi.fr/software/seaview).

## Plotting the phylogeny and the counts of *somA*-like loci

First prepare a folder for your results if it's not there already:

	$ mkdir results

To plot the results of the `2_FindSomAs` pipeline, I made a separate R script called `scripts/somACounts.R`. I ran this script in Rstudio, using the following conda environment (see the README at `2_FindSomAs` for more details on conda):

	$ mamba env create -f envs/r-ggtree-env.yml
	$ mamba activate r-ggtree

Then I opened Rstudio from within the environment:

	$ open -a RStudio

And obtained a figure in `results/tree_heatmap_v1.png`. As input it takes the `../2_FindSomAs/results/Loci_counts.txt` file. It also takes the newick tree produced above and a table mapping species to strains (`Sample_Species_map.txt`), both available in the directory `data` here.
