# load packages
library(treeio)    
library(ggtree)    
library(ape)       
library(phytools)  
library(phangorn)  
library(ggplot2)

# read trees
raxml_td <- read.tree("data/alignment_results/MAFFT_aligned.fasta.treefile")
iq_td <- read.tree("data/alignment_results/MAFFT_aligned.fasta.contree")

# check root location
plot(raxml_td)
nodelabels()
raxml_rooted = root(raxml_td, node=17, resolve.root=TRUE)
plot(raxml_rooted)

plot(iq_td)
nodelabels()
iq_rooted = root(iq_td, node=17, resolve.root = TRUE)
plot(iq_rooted)
