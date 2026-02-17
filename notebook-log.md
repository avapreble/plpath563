# Reproducible Notebook Log

## Question and Data

Regions with different first-line HIV treatments, and variations in the emergence of drug resistance

Has treatment resistance emerged multiple times independently in response to similar selective pressures, or does it exist more so as a result of transmission across regions?

I will use publicly available data of HIV-1 pol sequences from 3 regions - North America, Southern Africa, and Asia - and compare phylogenies. I selected these regions because they each have a different first-line antiretroviral therapy methods.  


## Quality Control

I do not have my data downloaded yet because I am still trying to figure out what specific files I should be downloading.
- update (2/16/2026): data downloaded, details below

I plan to use methods taught in class to assemble and align my data with as much accuracy as possible. 
- update (2/16/2026): going to run MAFFT on combined file of all HIV-1 sequences (5 random BLAST from each of the 3 chosen countries, combined_HIV2018)

I will maintain detailed notes throughout the process to ensure my work is reproducible, as well as to document steps taken in QC.


## Plan
1) download necessary/appropriate data from public database (.fasta files)
2) clean datasets/quality control
3) analyze alignments (MAFFT, MUSCLE?) and compare variation between each region's data
4) build tree and make inferences about how treatment resistance is most commonly spread in HIV-1 

## Downloading data
I downloaded my data from the Los Alamos National Laboratory HIV sequence database (https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html).
I input the following parameters in to the Sequence Search Interface:
- Sampling year: 2018
- Virus: HIV-1
- Subtype: any subtype
- Genomic region: Pol CDS
- Geographic region: *North America/Africa/Asia*

I then selected 5 random BLAST sequences from each of my 3 chosen countries - United States, Mozambique, and Thailand - of the most prominent subtype in each region.
### Accession numbers
- US (subtype B): MK169417, MK169418, MK169422, MK169441, MK169431
- Mozambique (subtype C): OK649266, OK649268, OK649270, OK649287, OK649295
- Thailand (subtype CRF01_AE): ON903098, ON903078, ON863057, ON863022, ON863108

## Aligning Data (one method for now, HW assignment) 

I chose to use MAFFT as my MSA method for this assignment


*MAFFT is a multiple sequence alignment program that works functions via 2 main techniques - (1) identification of homologous regions with Fast Fourier transform, (2) an efficient, quicker scoring system.*

*Assumptions*
- The chosen sequence regions are homologous and closely related
- Needleman-Wunsch algorithm for computing alignment without rearrangement

*Limitations*
- CPU times proportional to sequence length for more conserved regions, but closer to Needleman-Wunsch for more distantly related sequences
- More reliable guide trees produce more reliable results, an inaccurate guide tree can mess with alignment
- Fast Fourier Transform algorithm not efficient for more distantly related sequences 


- I made a combined file of all the downloaded sequences called *combined_HIV2018*

mafft combined_HIV2018.fasta > MAFFT_aligned.fasta   # code ran for running MAFFT (already installed) inside data folder of project

- This made the file *MAFFT_aligned.fasta* in the data folder of my project, which I put in to a sub folder of alignment results called *alignment_results*

