# running RAxML on alignment generated with MAFFT (15 taxa) 

#loading packages

library(ape)


# loading RAxML tree
tree <- read.tree("results/mafft_support.raxml.support")

#plot tree
plot(tree)
