# Reproducible Notebook Log

## Question and Data

Regions with different first-line HIV treatments, and variations in the emergence of drug resistance

Has treatment resistance emerged multiple times independently in response to similar selective pressures, or does it exist more so as a result of transmission across regions?

I will use publicly available data of HIV-1 pol sequences from 3 regions - North America, Southern Africa, and Asia - and compare phylogenies. I selected these regions because they each have a different first-line antiretroviral therapy methods.  


## Quality Control

I do not have my data downloaded yet because I am still trying to figure out what specific files I should be downloading.

I plan to use methods taught in class to assemble and align my data with as much accuracy as possible. 

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

I then selected 5 random BLAST sequences from each of my 3 chosen countries - United States, Monzambique, and Thailand
