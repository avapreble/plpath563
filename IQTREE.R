#loading packages

library(ape)

# load tree
iqtree <-read.tree('data/alignment_results/MAFFT_aligned.fasta.treefile')

# ladderize tree for visual clarity
iqtree_ladderize <- ladderize(iqtree)
plot(iqtree_ladderize, cex=0.5)
nodelabels(cex=0.7)


# I realigned my data for quality control purposes, removing columns with only gaps to see the impact
```
iqtree3 -s data/alignment_results/MAFFT_aligned_trimmed.fasta -m MFP -B 1000
```


#making tree with improved alignment
new_iqtree <- read.tree("data/alignment_results/MAFFT_aligned_trimmed.fasta.treefile")
#ladderize
new_iqtree_ladderize <- ladderize(new_iqtree)
plot(new_iqtree_ladderize, cex=0.5)
nodelabels(cex=0.7)





