# removing gap-heavy regions in alignment
library(ape)

# removing columns with only gaps
trim_alignment <- del.colgapsonly(dna)

# save improved alignment
write.dna(trim_alignment, file="data/alignment_results/MAFFT_aligned_trimmed.fasta", format="fasta")

