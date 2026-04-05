---
editor_options: 
  markdown: 
    wrap: 72
---

# **Reproducible Notebook Log**

## **Question and Data**

Regions with different first-line HIV treatments, and variations in the
emergence of drug resistance

Do prominent regional HIV-1 subtypes exhibit convergent evolution of
drug resistance mutations in the pol gene despite variations in
first-line antiretrovirals?

I will use publicly available data of HIV-1 pol sequences from 3
countries - United States, Mozambique, and Thailand - and compare
phylogenies. I selected these regions to ensure diverse data and because
of their variation in use of antiretroviral treatments.

## **Quality Control**

Using methods taught in class, I assembled and aligned my data with as much accuracy as possible. I ran MAFFT on my combined file (combined_HIV2018) of all HIV-1 pol sequences (5 random BLAST results from each of the 3 different countries).

I will maintain detailed notes throughout the process to ensure my work
is reproducible, as well as to document steps taken in QC.

## **Plan**

1)  download necessary/appropriate data from public database (.fasta
    files)
2)  clean datasets/quality control
3)  analyze alignments and compare variation between
    each region's data
4)  build tree and make inferences about how treatment resistance is
    most commonly spread in HIV-1

## **Downloading data**

I downloaded my data from the Los Alamos National Laboratory HIV
sequence database
(<https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html>). I input
the following parameters in to the Sequence Search Interface:

- Sampling year: 2018 
- Virus: HIV-1 
- Subtype: any subtype 
- Genomic region: Pol CDS 
- Geographic region: *North America/Africa/Asia*

I then selected 5 random BLAST sequences from each of my 3 chosen
countries - United States, Mozambique, and Thailand - of the most
prominent subtype in each region. 

**Accession numbers** 

- **US (subtype B):** MK169417, MK169418, MK169422, MK169441, MK169431 

- **Mozambique (subtype C):** OK649266, OK649268, OK649270, OK649287, OK649295 

- **Thailand (subtype CRF01_AE):** ON903098, ON903078, ON863057, ON863022, ON863108


## **Aligning Data MSA** 
#### (MAFFT, HW assignment)

### Description:
*MAFFT is a multiple sequence alignment program that works functions via
2 main techniques - (1) identification of homologous regions with Fast
Fourier transform, (2) an efficient, quicker scoring system.*

**Assumptions** - The chosen sequence regions are homologous and closely
related - Needleman-Wunsch algorithm for computing alignment without
rearrangement

**Limitations** - CPU times proportional to sequence length for more
conserved regions, but closer to Needleman-Wunsch for more distantly
related sequences - More reliable guide trees produce more reliable
results, an inaccurate guide tree can mess with alignment - Fast Fourier
Transform algorithm not efficient for more distantly related sequences

-   I made a combined file of all the downloaded sequences called
    *combined_HIV2018*

``` 
mafft combined_HIV2018.fasta > MAFFT_aligned.fasta   # code ran for running MAFFT (already installed) inside data folder of project
```

-   This made the file *MAFFT_aligned.fasta* in the data folder of my
    project, which I put in to a sub folder of alignment results called
    *alignment_results*
    
    
    
## **Building simpler trees to start**

#### loaded packages:
``` 
library(ape)
library(adegenet)
library(phangorn)
```

### Distance-based tree
**Distance-based trees are built via computation of pairwise genetic distances and hierarchical clustering algorithms, in this case Neighbor-Joining.**

**Assumptions:** through choosing the TN93 model, we assume there are varying rates of transitions and transversions and heterogeneous base frequencies.

**Limitations:** only outputs 1 tree when there could be multiple equally-possible options. Also, the optimization of pairwise genetic distances and the built tree is very algorithmic rather than explicitly based in evolutionary biology. 

#### Loading sample data

``` 
dna <- fasta2DNAbin(file="data/alignment_results/MAFFT_aligned.fasta")
``` 

#### Computation of genetic distances (Tamura and Nei 1993 model)

*(idk if this is correct model to use it was this that we used in class demonstration though, may edit)*


``` 
D <- dist.dna(dna, model="TN93")
```


#### getting NJ tree

``` 
tre <- nj(D)
```


#### ladderize

``` 
tre <- ladderize(tre)
```


#### plot tree
``` 
plot(tre, cex=.6)
title("Simple NJ tree")
```


### Parsimony-based tree

**Parsinomy-based trees are built upon the basis of using the smallest number of genetic sequence changes between each taxa.**

**Assumptions:** the simplest scenario, with the least amount of genetic changes, is the most likely. This is more accurate when the genetic changes are less compared to something rapidly evolving and mutating

**Limitations:** assumes the simplest scenario, which doesn't account for large eveolutionary events well. This model is also more simple and computer intensive than some of it's counterparts.   

#### R packages loaded above

#### Conversion of sample data in to phangorn object
``` 
dna2 <- as.phyDat(dna)
```

#### Starting tree 

``` 
starting_tree <- nj(dist.dna(dna,model="raw"))
parsimony(starting_tree, dna2)
```

#### Searching for tree with max. parsimony

``` 
parsimony_tree <- optim.parsimony(starting_tree, dna2)
```

#### 
- **console: final p-score 1037 after 0 nni operations**


#### Root the tree:
``` 
rtre = root(parsimony_tree, node = 19)
```
#### Tree plot

``` 
plot(parsimony_tree, cex=0.6)
title("Parsimony tree")
```




# **Maximum Likelihood trees and inference**

## **RAxML-NG v1.2.2**


#### Input alignment:

data/alignment_results/MAFFT_aligned.fasta


#### Verify format command:
```
raxml-ng --check --msa data/alignment_results/MAFFT_aligned.fasta --model GTR+G --prefix results/mafft_check
```
*RAxML detected two identical sequences, so it reduced the alignment for computational purposes (C.MZ.2018.S04.OK649268 and C.MZ.2018.S06.OK649270)*
*Alignment has 3094 sites and 15 taxa*
 
 

#### Reduced alignment for analysis

results/mafft_check.raxml.reduced.phy


#### Build RAxML ML tree command

```
raxml-ng --search --msa results/mafft_check.raxml.reduced.phy --model GTR+G --prefix results/mafft_ml
```
*tree saves as results/mafft_ml.raxml.bestTree*

#### Output: 

- logLikelihood: -9274.482730

- 14 taxa due to duplicate removed for computation

# 


#### Bootstrap trees command:
```
raxml-ng --bootstrap --msa results/mafft_check.raxml.reduced.phy --model GTR+G --prefix results/mafft_boot
```


#### Bootstrap replicate tree output:

results/mafft_boot.raxml.bootstraps



#### Map bootstrap to best ML tree command:
```
raxml-ng --support --tree results/mafft_ml.raxml.bestTree --bs-trees results/mafft_boot.raxml.bootstraps --prefix results/mafft_support
```
#### Tree:
results/mafft_support.raxml.support


#### Commands: 

loading packages:
```
library(ape)
```

loading RAxML tree:
```
tree <- read.tree("results/mafft_support.raxml.support")
```

plot tree:
``` 
plot(tree)
```

# 
# 

## **IQ-TREE v3.0.1**

### Description:
IQ-TREE is a software used for phylogenetic inference by maximum likelihood.The software package includes model selection through ModelFinder, an effective search algorithm, and fast bootstrapping. 

**Assumptions:**
(https://iqtree.github.io/doc/Assessing-Phylogenetic-Assumptions) 

- treelikeness: all sites in aligned data were yielded from the same tree

- stationarity: constant frequency of nucleotides and amino acids through time

- reversibility: substitutions equally likely in both directions

- homogeneity: constant substitution rate through time


**Limitations:**

- quality of alignment determines results

- extensive bootstrapping still presents uncertainty

- weak phylogenetic signals within a dataset may require parameter adjustment 


### Input alignment:

data/alignment_results/MAFFT_aligned.fasta 

# 

### Command input
```
~/Documents/software/iqtree-3.0.1-macOS/bin/iqtree3 -s data/alignment_results/MAFFT_aligned.fasta -m MFP -B 1000
```

### Output 
**MAFFT_aligned.fasta.log**

- 15 sequences, 3094 sites

- 1 pair of identical seq (NOTE: C.MZ.2018.S06.OK649270 is identical to C.MZ.2018.S04.OK649268 but kept for subsequent analysis
Checking for duplicate sequences: done in 6.29425e-05 secs using 66.73% CPU)

- best-fit model (Bayesian Information Criterion): TIM3+F+I



## MrBayes - Bayesian Inference

### Description: xxxxxxxxxxx

**Assumptions**
- "Concatenation methods implicitly assume that all gene loci share the same topology and branch lengths." (Huelsenbeck and Ronquist)
- data homogeneity 

**Limitations**
- dependent on correctness of model


### Commands

#### First, converting .fasta to .nex
```
library(ape)
dna <- read.dna("data/alignment_results/MAFFT_aligned.fasta", format = "fasta")
write.nexus.data(dna, file = "data/alignment_results/MAFFT_aligned.nex", format= "dna")
```

#### Running MrBayes
```
mb

MrBayes > execute data/alignment_results/MAFFT_aligned.nex

MrBayes > lset nst=6 rates=propinv

mcmc ngen=10000

sump 

sumt 
```

#### Output

- summarization of all sample trees to one summary tree

- confidence values (posterior probabilities)


