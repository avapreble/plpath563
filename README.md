# **Phylogenetic Analysis of HIV-1** ***pol*** **Sequences From 3 Geographical Regions**


#### Ava Preble | PLPATH563: Phylogenetic Analysis of Molecular Data 

---


***QUICK LINKS TO FINAL PROJECT:***

- [FINAL RESEARCH PAPER](FinalResearchPaper_plpath563.pdf)
- [REPRODUCIBILITY SCRIPT](notebook-log.md)


---

#### This project uses phylogenetic methods (including distance and parsimony, maximum likelihood, and Bayesian inference) to analyze HIV-1 *pol* sequences from three regions: the United States, Mozambique, and Thailand. 


Data was gathered from Los Alamos National Laboratory HIV sequence database. 
Steps for reproducibility can be found in the [notebook-log markdown file](notebook-log.md). 

Required software: *R, RStudio, MAFFT, IQ-TREE, RAxML-NG, MrBayes, BEAST*

-----

#### **REPOSITORY NAVIGATION**

##### **MAIN DIRECTORY: PLPATH563** 
***plpath563/***  --> *README and reproducibility notebook files, final research paper, R project, project folders (SEE BELOW)*

- **data/** 
individual FASTA files by country, combined FASTA file, folder with MAFFT outputs

  - **data/alignment_results/** 
    - MAFFT outputs  
  


- **figures/**  
Constructed trees and Tracer graphs (png or pdf)

- **scripts/**  
individual R scripts for each phylogenetic method used

- **results/**  
RAxML outputs, folder with distance matrix

  - **results/matrix/**  
    - distance matrix CSV


- **beast/**  
BEAST and BEAUti2 analyses

- **other/**  
unrelated homework and in-class practices