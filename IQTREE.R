#loading packages

library(ape)

# load tree
iqtree <-read.tree('data/alignment_results/MAFFT_aligned.fasta.treefile')
plot(iqtree)
