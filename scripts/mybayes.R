# conversion .fasta to .nex
library(ape)
dna <- read.dna("data/alignment_results/MAFFT_aligned.fasta", format = "fasta")
write.nexus.data(dna, file = "data/alignment_results/MAFFT_aligned.nex", format= "dna")

# MrBayes

mb

MrBayes > execute data/alignment_results/MAFFT_aligned.nex

MrBayes > lset nst=6 rates=propinv

mcmc ngen=10000

sump (summarizes parameter estimates)

sumt (summarizes trees)