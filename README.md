# Convergence event counting and probability calculation (conv_cal)

Author: Zhengting Zou

Version 0.3 (2021-03-13)

__Description__: This is an analysis pipeline for counting the observed number of 
convergence events within gene sequence alignments between N > 1 branches 
on a phylogeny, and calculating the expected total probability of parallel and 
convergent events accordingly. The method is used in:

> Zou Z, and Zhang J. 2015. Are convergent and parallel amino acid substitutions 
in protein evolution more prevalent than neutral expectations? 
Mol. Biol. Evol., 32: 2085-2096.

This is an updated version of the `convCal` package previously distributed, and 
instead of analyzing convergence between two branches, it can now derive the 
observed number and expected probability of convergence between N (N >= 2) 
branches.

### Required softwares and packages:
 - PAML (version >= 4.9j)
 - Python (version >= 3.7)
 - Python packages: numpy, scipy, ete3

### Required input files:
(here stored in `00_input/`)
- Sequence alignment file of fasta format that PAML can accept.
- A .dat file defining the amino acid substitution model used by PAML.
- A newick style tree file describing the phylogeny of all species concerned.

### Analysis step 1: Inference of ancestral sequences and parameters by PAML
(here stored in `01_ancestral/`)

__[script usage]__

`./run_paml.py <input directory> <output directory> <path to tree file> `
`<path to .dat file> [<number of CPU cores to use>, ]`

Command-line example stored in `cmd.bash`.

All .fasta and .fa files in the input directory will be subject to a 
`codeml `analysis, the result of which stored in a sub-directory in the 
output directory. A tree file named `tree_labeled.nw` will be generated 
for the convenience of designating focal branch groups in the next step.

### Analysis step 2: Counting and calculating convergence
(here stored in `02_convergence/`)

__[script usage]__

`./calc_expconv.py <input directory> <output directory> <path to group file> `
`<path to .dat file> [<a.a. frequency mode>, ]`

Command-line example stored in `cmd.bash`.

Use the directory that contain `codeml` results as input directory. 

The group file contain one or more lines, each designating a group of 
branches between which convergence is to be counted and calculated, 
separated by whitespaces. Use the end point node name as the name 
of a branch, according the `tree_labeled.nw` generated in the previous 
step (note that the number prefix, e.g. "7\_"  is not needed for a leaf 
node). Example of a group file: `branch_group.list`.

The last command-line parameter can be `site` or `gene`, corresponding 
to the JTT-f(site) and JTT-f(gene) model in the original study, or it can be 
a frequency vector / matrix file that can be read by `numpy.loadtxt`. 

### Output files:

Each gene (with a single `.fasta` file as input into the pipeline) will 
generate an `.obsconv.tsv` file and an `.expconv.tsv` file.

The `.obsconv.tsv` contain multiple lines corresponding to the 
branch groups analyzed. Each line contain the branch group (branches 
separated by ","), the number of parallel and convergent events 
counted across the gene, and detail of the events in parentheses 
(site index and ancestral / derived amino acid states of the branches). 


Similarly, each line of `.expconv.tsv` contain the branch group, 
probability sums of paralle / convergent event across the whole gene 
(i.e. expected numbers), and site-specific probabilities.

