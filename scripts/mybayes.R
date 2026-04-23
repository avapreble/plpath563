# conversion .fasta to .nex (pre-realignment)
library(ape)
dna <- read.dna("data/alignment_results/MAFFT_aligned.fasta", format = "fasta")
write.nexus.data(dna, file = "data/alignment_results/MAFFT_aligned.nex", format= "dna")

# MrBayes

# Running MrBayes in terminal, begin inside project folder (plpath563)

````
mb

MrBayes > execute data/alignment_results/MAFFT_aligned.nex

MrBayes > lset nst=6 rates=propinv

mcmc ngen=500000 

sump #summarizes parameter estimates

sumt #summarizes trees
````

# visualize
mrbayes_tree <- read.nexus("data/alignment_results/MAFFT_aligned.nex.con.tre") #loading in MB tree 
plot(mrbayes_tree, cex=0.7)
nodelabels(cex=0.8) #posterior probabilities




# MrBayes with trimmed alignment

# conversion from .fasta to .nex
mrbayes_trimmed <- read.dna("data/alignment_results/MAFFT_aligned_trimmed.fasta", format = "fasta")
write.nexus.data(mrbayes_trimmed, file = "data/alignment_results/MAFFT_aligned_trimmed.nex", format= "dna")




# running MrBayes with trimmed alignment (in terminal)
```
mb

MrBayes > execute data/alignment_results/MAFFT_aligned_trimmed.nex

mcmc ngen=500000 

sump #summarizes parameter estimates

sumt #summarizes trees
```

#visualize
mb_trimmed_tree <- read.nexus("data/alignment_results/MAFFT_aligned_trimmed.nex.con.tre")
mb_trimmed_ladderized <- ladderize(mb_trimmed_tree)
plot(mb_trimmed_ladderized, cex=0.7)
nodelabels(cex=0.8)






