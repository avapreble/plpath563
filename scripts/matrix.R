library(ape) 

dna <- read.dna("data/alignment_results/MAFFT_aligned.fasta", format="fasta")
distmatrix <- dist.dna(dna, model="raw")
as.matrix(distmatrix)


# saving matrix as .csv spreadsheet
write.csv(as.matrix(distmatrix), file = "results/matrix/distmatrix.csv")

