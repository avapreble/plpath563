### Distance and parsimony methods on data (HW)

# loading R packages
library(ape)
library(adegenet)
library(phangorn)

## Distance-based tree

# Distance-based trees are built via computation of pairwise genetic distances and hierarchical clustering algorithms, in this case Neighbor-Joining.

# Assumptions: through choosing the TN93 model, we assume there are varying rates of transitions and transversions and heterogeneous base frequencies.

# Limitations: only outputs 1 tree when there could be multiple equally-possible options. Also, the optimization of pairwise genetic distances and the built tree is very algorithmic rather than explicitly based in evolutionary biology. 

# loading sample data
dna <- fasta2DNAbin(file="data/alignment_results/MAFFT_aligned.fasta")

# computation of genetic distances (Tamura and Nei 1993 model)
# (idk if this is correct model to use it was this that we used in class demonstration though, may edit)
D <- dist.dna(dna, model="TN93")

# getting NJ tree
tre <- nj(D)

#ladderize
tre <- ladderize(tre)

#plot tree
plot(tre, cex=.6)
title("Simple NJ tree")



## Parsimony-based tree

# Parsinomy-based trees are built upon the basis of using the smallest number of genetic sequence changes between each taxa.

# Assumptions: the simplest scenario, with the least amount of genetic changes, is the most likely. This is more accurate when the genetic changes are less compared to something rapidly evolving and mutating
# Limitations: assumes the simplest scenario, which doesn't account for large eveolutionary events well. This model is also more simple and computer intensive than some of it's counterparts.   


# R packages loaded above

# conversion of sample data in to phangorn object
dna2 <- as.phyDat(dna)

# starting tree 
starting_tree <- nj(dist.dna(dna,model="raw"))
parsimony(starting_tree, dna2)

# searching for tree with max. parsimony
parsimony_tree <- optim.parsimony(starting_tree, dna2)
# console: final p-score 1037 after 0 nni operations(??)

# root the tree:
rtre = root(parsimony_tree, node = 19)

# tree plot
plot(rtre, cex=0.6)
title("Parsimony tree")



