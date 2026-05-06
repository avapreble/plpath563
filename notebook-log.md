---
editor_options: 
  markdown: 
    wrap: 72
---

# **Reproducible Notebook Log**

---

#### Software used:
- RStudio (Version 2025.05.1+513)
- R (Version 4.5.1)
- MAFFT (Version 7.526)
- IQ-TREE (Version 3.0.1)
- RAxML-NG (Version 1.2.2)
- MrBayes (Version 3.2.7a)
- BEAST (Version 2.7.7)

---

## **Question and Data**

### ***How do phylogenetic clustering patterns of HIV-1 pol sequences differ across three distinct geographic regions, and what do these patterns suggest about viral evolution in regions with varying treatment contexts?***

#### Looking at phylogenetic clustering of HIV-1 *pol* sequences from three geographic regions.


I will use publicly available data of HIV-1 pol sequences from 3
countries - United States, Mozambique, and Thailand - and compare
phylogenies. I selected these regions to ensure diverse data and because
of their variation in use of antiretroviral treatments.

------------------------------------------------------------------------

### **Step 1: Obtaining Data**

The raw data used in this project was obtained from the Los Alamos
National Laboratory HIV sequence database
(<https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html>).

Sequence Search Interface parameters input:

-   Sampling year: 2018
-   Virus: HIV-1
-   Subtype: any subtype
-   Genomic region: Pol CDS
-   Geographic region: *North America/Africa/Asia*

I then selected 5 random BLAST sequences from each of my 3 chosen
countries - United States, Mozambique, and Thailand - of the most
prominent subtype in each region.

**Accession numbers**

-   **US (subtype B):** MK169417, MK169418, MK169422, MK169441, MK169431

-   **Mozambique (subtype C):** OK649266, OK649268, OK649270, OK649287,
    OK649295

-   **Thailand (subtype CRF01_AE):** ON903098, ON903078, ON863057,
    ON863022, ON863108

------------------------------------------------------------------------

### **Step 2: MAFFT for Multiple Sequence Alignment**


*MAFFT is a multiple sequence alignment program that functions via
2 main techniques - (1) identification of homologous regions with Fast
Fourier transform, (2) an efficient, quicker scoring system.*

**Assumptions**

- The chosen sequence regions are homologous and closely
related 

- Needleman-Wunsch algorithm for computing alignment without
rearrangement

**Limitations** 

- CPU times proportional to sequence length for more conserved regions, but closer to Needleman-Wunsch for more distantly related sequences 

- More reliable guide trees produce more reliable results, an inaccurate guide tree can mess with alignment 

- Fast Fourier Transform algorithm not efficient for more distantly related sequences

3 .fasta files with each country's 5 pol sequences were combined in to a
collective .fasta before MAFFT alignment (*combined_HIV2018.fasta*).

Running MAFFT in terminal to create alignment file:

```         
mafft combined_HIV2018.fasta > MAFFT_aligned.fasta   
```

The output was *MAFFT_aligned.fasta*, which was put in:

```         
data/alignment_results
```

.

***Quality control note:*** *Using methods taught in class, I assembled
and aligned my data with as much accuracy as possible. I ran MAFFT on my
combined file (combined_HIV2018) of all HIV-1 pol sequences (5 random
BLAST results from each of the 3 different countries). Later, I decided
to try to trim the data to see if it had any significant effect on the
analyses. After assessing IQ-TREE and MrBayes outputs with the trimmed
alignment and comparing them to the originals, I concluded I would not
rework everything with the re-aligned data, since difference in results
was minuscule. See the creation of the trimmed alignment to follow:*

##### Making trimmed alignment (quality control)

[R scripts here](script/alignmentQC.R) I created my original alignment
and analyzed without trimming, but then I made a second alignment of the
same data that was trimmed, removing columns that only consisted of gaps
(see below). This did not change overall topology or support values in a
major way.

```r
library(ape)

# removing columns with only gaps
trim_alignment <- del.colgapsonly(dna)

# save improved alignment
write.dna(trim_alignment, file="data/alignment_results/MAFFT_aligned_trimmed.fasta", format="fasta")
```

-   Similar to the original, we now have *MAFFT_aligned_trimmed.fasta*
    in *alignment_results*

------------------------------------------------------------------------

### **Step 3: Distance & Parsimony Trees**

[R scripts here](scripts/DistanceAndParsimony.R)

#### Loaded R packages:

```r
library(ape)
library(adegenet)
library(phangorn)
```

### Distance-based tree

**Distance-based trees are built via computation of pairwise genetic
distances and hierarchical clustering algorithms, in this case
Neighbor-Joining.**

**Assumptions:** through choosing the TN93 model, we assume there are
varying rates of transitions and transversions and heterogeneous base
frequencies.

**Limitations:** only outputs 1 tree when there could be multiple
equally-possible options. Also, the optimization of pairwise genetic
distances and the built tree is very algorithmic rather than explicitly
based in evolutionary biology.

```r
# R packages loaded above


# Loading sample data


dna <- fasta2DNAbin(file="data/alignment_results/MAFFT_aligned.fasta")

# Computation of genetic distances (Tamura and Nei 1993 model)

# (I used this model because it is what we used in class)


D <- dist.dna(dna, model="TN93")


# Getting NJ tree


tre <- nj(D)


# Ladderize


tre <- ladderize(tre)


# Plot tree


plot(tre, cex=.6)
title("Simple NJ Distance tree")
```

### Parsimony-based tree

**Parsimony-based trees are built upon the basis of using the smallest
number of genetic sequence changes between each taxa.**

**Assumptions:** the simplest scenario, with the least amount of genetic
changes, is the most likely. This is more accurate when the genetic
changes are less compared to something rapidly evolving and mutating

**Limitations:** assumes the simplest scenario, which doesn't account
for large evolutionary events well. This model is also more simple and
computer intensive than some of it's counterparts.

```r
#### R packages loaded above

# Conversion of sample data in to phangorn object


dna2 <- as.phyDat(dna)


# Starting tree


starting_tree <- nj(dist.dna(dna,model="raw"))
parsimony(starting_tree, dna2)


# Searching for tree with max. parsimony


parsimony_tree <- optim.parsimony(starting_tree, dna2)

# console: *final p-score 1037 after 0 nni operations*


# Root the tree:


rtre = root(parsimony_tree, node = 19)

# Tree plot


plot(parsimony_tree, cex=0.6)
title("Parsimony tree")
```

------------------------------------------------------------------------

### **Step 4: Maximum Likelihood (RAxML, IQ-TREE)**

### **RAxML-NG v1.2.2**

[R scripts here](scripts/raxmlHIV.R)

#### Input alignment:

```         
data/alignment_results/MAFFT_aligned.fasta
```

#### Verify format command:

```         
raxml-ng --check --msa data/alignment_results/MAFFT_aligned.fasta --model GTR+G --prefix results/mafft_check
```

*RAxML detected two identical sequences, so it reduced the alignment for
computational purposes (C.MZ.2018.S04.OK649268 and
C.MZ.2018.S06.OK649270)* 

*Alignment has 3094 sites and 15 taxa*

#### Reduced alignment for analysis

```         
results/mafft_check.raxml.reduced.phy
```

#### Build RAxML ML tree command

```         
raxml-ng --search --msa results/mafft_check.raxml.reduced.phy --model GTR+G --prefix results/mafft_ml
```

*tree saves as* `results/mafft_ml.raxml.bestTree`

#### Output:

-   logLikelihood: -9274.482730

-   14 taxa due to duplicate removed for computation

# 

#### Bootstrap trees command:

```         
raxml-ng --bootstrap --msa results/mafft_check.raxml.reduced.phy --model GTR+G --prefix results/mafft_boot
```

#### Bootstrap replicate tree output: `results/mafft_boot.raxml.bootstraps`

#### Map bootstrap to best ML tree command:

```         
raxml-ng --support --tree results/mafft_ml.raxml.bestTree --bs-trees results/mafft_boot.raxml.bootstraps --prefix results/mafft_support
```

#### Tree: `results/mafft_support.raxml.support`

#### Commands:

```r
# Loading packages:


library(ape)


# Loading RAxML tree:


tree <- read.tree("results/mafft_support.raxml.support")


# Plot tree:


plot(tree)
```


#### Visualize tree with ggtree:

```r
#loading packages
library(ggtree)    
library(ape)         
library(ggplot2)
library(treeio)


# loading RAxML tree
raxmltree <- read.tree("results/mafft_support.raxml.support")
raxmltree <- ladderize(raxmltree)


# tip labels, colors, formating
colorRegionsraxml <- data.frame(label = raxmltree$tip.label)
colorRegionsraxml$Region <- ifelse(grepl("C.MZ", colorRegionsraxml$label), "Mozambique", 
                              ifelse(grepl("TH", colorRegionsraxml$label), "Thailand", "USA"))

colorRegionsraxml$new_label <- ave(colorRegionsraxml$Region, colorRegionsraxml$Region,
                              FUN = function(x) paste0(x, "_", seq_along(x)))

#plot tree
raxmltreeplot <- ggtree(raxmltree) %<+% colorRegionsraxml + 
  geom_tiplab(aes(label = new_label, color = Region), size = 3) +
  xlim(0, max(ggtree(raxmltree) $data$x) * 1.4) +
  geom_text2(aes(label = label, subset = !isTip), 
           size = 2.5, hjust = -2.5) 


raxmltreex <- raxmltreeplot + 
  theme_tree2() + 
  theme(legend.position = "left", text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(label = "•", size = 10))) 

raxmltreex

# save image
ggsave("figures/raxmltree1.png", plot = raxmltreex, width = 8, height = 6, dpi = 300)
```



## **IQ-TREE v3.0.1**

[R scipt here](scripts/IQTREE.R)

### Description:

IQ-TREE is a software used for phylogenetic inference by maximum
likelihood.The software package includes model selection through
ModelFinder, an effective search algorithm, and fast bootstrapping.

**Assumptions:**
(<https://iqtree.github.io/doc/Assessing-Phylogenetic-Assumptions>)

-   treelikeness: all sites in aligned data were yielded from the same
    tree

-   stationarity: constant frequency of nucleotides and amino acids
    through time

-   reversibility: substitutions equally likely in both directions

-   homogeneity: constant substitution rate through time

**Limitations:**

-   quality of alignment determines results

-   extensive bootstrapping still presents uncertainty

-   weak phylogenetic signals within a dataset may require parameter
    adjustment

#### Input alignment:

```         
data/alignment_results/MAFFT_aligned.fasta 
```

#### Command input:

```         
~/Documents/software/iqtree-3.0.1-macOS/bin/iqtree3 -s data/alignment_results/MAFFT_aligned.fasta -m MFP -B 1000
```

#### Output:

**MAFFT_aligned.fasta.log**

-   15 sequences, 3094 sites

-   1 pair of identical seq (NOTE: C.MZ.2018.S06.OK649270 is identical
    to C.MZ.2018.S04.OK649268 but kept for subsequent analysis Checking
    for duplicate sequences: done in 6.29425e-05 secs using 66.73% CPU)

-   best-fit model (Bayesian Information Criterion): TIM3+F+I

#### Visualize with ggtree:

```r
library(ggtree)    
library(ape)         
library(ggplot2)
library(treeio)

# load tree
iqtree <- read.tree("data/alignment_results/MAFFT_aligned.fasta.treefile") 


# tip labels, colors, formating
colorRegions <- data.frame(label = iqtree$tip.label)
colorRegions$Region <- ifelse(grepl("C.MZ", colorRegions$label), "Mozambique", 
                              ifelse(grepl("TH", colorRegions$label), "Thailand", "USA"))

colorRegions$new_label <- ave(colorRegions$Region, colorRegions$Region,
                              FUN = function(x) paste0(x, "_", seq_along(x)))

iqtreePlot <- ggtree(iqtree) %<+% colorRegions +
  geom_tiplab(aes(label = new_label, color = Region), size = 3) +
  geom_nodelab(aes(label = label), size = 2, hjust = -0.5, vjust = -0.3) +
  xlim(0, max(ggtree(iqtree)$data$x) * 1.4)


iqtreePlotx <- iqtreePlot + theme_tree2() + 
  theme(legend.position = "left", text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(label = "•", size = 10))) 

iqtreePlotx

# save image
ggsave("figures/iqtree1.png", plot = iqtreePlotx, width = 8, height = 6, dpi = 300)
```

### Re-running IQ-TREE on trimmed alignment for quality control

#### Command input:

```         
iqtree3 -s data/alignment_results/MAFFT_aligned_trimmed.fasta -m MFP -B 1000
```

### Output:

**MAFFT_aligned_trimmed.fasta.log**

#### Visualize:

```r
# making tree with improved alignment
new_iqtree <- read.tree("data/alignment_results/MAFFT_aligned_trimmed.fasta.treefile") 

#ladderize
new_iqtree_ladderize <- ladderize(new_iqtree)

plot(new_iqtree_ladderize, cex=0.5)
nodelabels(cex=0.7)
```
### **Step 5: Bayesian Inference (MrBayes, BEAST)**

## MrBayes 

[R scripts here](mrbayes.R)

### Description:

MrBayes is a software used for Bayesian inference modeling. Bayesian
inference begins with a prior probability that is continually updated
with the addition of more data to produce a posterior probability.

**Assumptions:**

-   "Concatenation methods implicitly assume that all gene loci share
    the same topology and branch lengths." (Huelsenbeck and Ronquist)

-   data homogeneity

-   correct substitution model chosen

**Limitations:**

-   dependent on alignment quality

-   dependent on correctness of model

#### Commands:

```r
# First, converting .fasta to .nex
library(ape)
dna <- read.dna("data/alignment_results/MAFFT_aligned.fasta", format = "fasta")
write.nexus.data(dna, file = "data/alignment_results/MAFFT_aligned.nex", format= "dna")
```

#### Running MrBayes:

```      
mb

MrBayes > execute data/alignment_results/MAFFT_aligned.nex

MrBayes > lset nst=6 rates=propinv

mcmc ngen=10000

sump 

sumt 
```

#### Output:

-   summarization of all sample trees to one summary tree

-   confidence values (posterior probabilities)

#### To visualize:

```r       
mrbayes_tree <- read.nexus("data/alignment_results/MAFFT_aligned.nex.con.tre") #loading in MB tree 
plot(mrbayes_tree, cex=0.7)
nodelabels(cex=0.8) #posterior probabilities
```

#### Visualizing with ggtree:
 
```r        
# load packages
library(ggtree)    
library(ape)         
library(ggplot2)
library(treeio)



# MrBayes - read in tree
MBtree <- read.nexus("data/alignment_results/MAFFT_aligned.nex.con.tre")
MBtreex <- read.mrbayes("data/alignment_results/MAFFT_aligned.nex.con.tre")

# base tree
plot(MBtree)
nodelabels()

# tip labels
colorRegions <- data.frame(label = MBtree$tip.label)
colorRegions$Region <- ifelse(grepl("C.MZ", colorRegions$label), "Mozambique", 
                              ifelse(grepl("TH", colorRegions$label), "Thailand", "USA"))

colorRegions$new_label <- ave(colorRegions$Region, colorRegions$Region,
                              FUN = function(x) paste0(x, "_", seq_along(x)))

# aesthetics

mbTreePlot <- ggtree(MBtreex) + geom_tiplab()
mbTreePlot <- ggtree(MBtreex) %<+% colorRegions + 
  
  geom_tiplab(aes(label = new_label, color = Region), size = 3) + 
  
  geom_text2(aes(label = sprintf("%.2f", as.numeric(prob)), subset = !isTip), 
             size = 2.5, hjust = -2.5) +
  
  xlim(0, max(ggtree(MBtreex)$data$x) * 1.3) + 
  
  
  theme_tree2() + 
  theme(legend.position = "left", text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(label = "•", size = 10))) 

mbTreePlot

# save image
ggsave("figures/mbtree1.png", plot = mbTreePlot, width = 8, height = 6, dpi = 300)
```

### Re-running MrBayes on trimmed alignment for quality control

#### Converting trimmed .fasta to .nex:

```         
mrbayes_trimmed <- read.dna("data/alignment_results/MAFFT_aligned_trimmed.fasta", format = "fasta")
write.nexus.data(mrbayes_trimmed, file = "data/alignment_results/MAFFT_aligned_trimmed.nex", format= "dna")
```

#### Running MrBayes on trimmed alignment:

```         
mb

MrBayes > execute data/alignment_results/MAFFT_aligned_trimmed.nex

mcmc ngen=500000 

sump #summarizes parameter estimates

sumt #summarizes trees
```

#### To visualize:

```         
mb_trimmed_tree <- read.nexus("data/alignment_results/MAFFT_aligned_trimmed.nex.con.tre")
mb_trimmed_ladderized <- ladderize(mb_trimmed_tree)
plot(mb_trimmed_ladderized, cex=0.7)
nodelabels(cex=0.8)
```


## BEAST and MCMC

BEAST is a software package that uses MCMC for Bayesian analysis of
molecular sequences. It estimates evolutionary parameters. The aligned
sequences are put in to BEAUti, which generates the .xml to input in
BEAST. BEAST outputs .log and .trees files. Then, the .log file
(beast/beast_BEAUti2_analysis.log) is loaded in to Tracer, which
presents visuals of posterior trace, likelihood, ESS, etc.

-   The figures created in Tracer (posterior and likelihood plots) after
    uploading the BEAUti2 output show MCMC convergence after the initial
    calibration period. The effective sample size (ESS) values, however,
    were low (posterior = 58, likelihood = 53), indicating that the
    various pol sequences were very similar to one another and the
    sequences weren't really long enough to generate a reliable result.

[posterior plot](figures/posterior_tracer_beast.pdf)

[likelihood plot](figures/likelihood_tracer_beast.pdf)

------------------------------------------------------------------------

### **Step 6: Distance Matrix**

#### Pairwise distances

```r
# Load R package needed:

         
library(ape) 


# Creating distance matrix:

dna <- read.dna("data/alignment_results/MAFFT_aligned.fasta", format="fasta")
distmatrix <- dist.dna(dna, model="raw")
as.matrix(distmatrix)


Saving matrix as .csv (spreadsheet)



        
write.csv(as.matrix(distmatrix), file = "results/matrix/distmatrix.csv")
```
[matrix here](results/matrix/distmatrix.csv)

-   The created distance matrix is comparing each of the 15 sequences to
    one another. Lower values indicate more similarity and higher values
    indicate less. Within each country, the strains have very high
    similarity. However, between each country, the strains are more
    variable.

------------------------------------------------------------------------

### Result interpretations

#### Plotting trees with ggtree

-   ggtree is an R package for visualizing phylogenetic trees and other
    relevant data in a more aesthetic, clean format than just using
    plot(tree).
