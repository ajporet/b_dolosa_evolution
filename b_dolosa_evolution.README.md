This repository contains the code and analyses used in:


### “De novo mutations mediate phenotypic switching in an opportunistic human lung pathogen ”<a id="de-novo-mutations-mediate-phenotypic-switching-in-an-opportunistic-human-lung-pathogen-"></a>

Alexandra J. Poret, Matthew Schaefers, Christina Merakou, Kathryn E. Mansour, Georgia K. Lagoudas, Ashley R. Cross, Joanna B. Goldberg, Roy Kishony, Ahmet Z. Uluer, Alexander J. McAdam, Paul C. Blainey, Sara O. Vargas, Tami D. Lieberman^, Gregory P. Priebe^ (^co-senior authors)


# Data description:<a id="data-description"></a>

In this article, we profile a three-person outbreak of _Burkholderia dolosa_, an opportunistic pathogen that infects those with cystic fibrosis. This cluster comprised 3 adults: Subject J (the suspected index patient), Subject Q (co-worker of Patient J), and Subject R (sibling of Patient Q). In this repo, Subject J is referred to as “P02”, Subject Q as “P06”, and Subject J as “P07.” Whole genome sequences from 987 isolates obtained from these patients are stored under BioProject number PRJNA1063312. Note that these sequences are published post-cutadapt and sickle processing (see paper methods), so the initial filtering steps may be skipped for those repeating these analyses. 

In this paper, we also reanalyze 112 isolates from Lieberman and Michel et al.’s 2011 work: Parallel bacterial evolution within multiple patients identifies candidate pathogenicity genes (doi: <https://doi.org/10.1038/ng.997>). Raw sequencing data for these isolates can be found under Bioproject number PRJNA255914. 


# Data Processing: <a id="data-processing"></a>

This repo is broken down into five subsections:

1. lab\_scripts\_and\_files

2. LPS\_phenotyping

3. main\_figures

4. sample\_processing

5. supplemental\_figures

6. supplemental\_tables


### Lab\_scripts\_and\_files<a id="lab_scripts_and_files"></a>

Lab\_scripts\_and\_files contains standardized scripts and functions utilized throughout the Lieberman lab. These functions are commonly called throughout this repo. We recommend downloading this folder and linking it to your PATH before running our analyses. 


## Sample\_processing<a id="sample_processing"></a>

Sample\_processing contains code used to process _B. dolosa_ WGS raw reads. This directory contains two subsections: JQR\_evolution (processing of newly collected isolates) and Lieberman\_et\_al\_2011\_reanalysis (processes Lieberman and Michel’s historical _B. dolosa_ isolates). Within these directories, “extract\_data\_from\_samples” details how raw reads are filtered (can be skipped using data downloaded from NCBI), “create\_candidate\_mutation\_table” tabulates candidate SNVs, and “create\_tree” assembles those SNVs into a phylogeny. “Breseq” contains all code needed to search for indels within _B. dolosa_.  


## Main\_figures <a id="main_figures"></a>

All scripts used to generate a figure are compiled here, along with their necessary data structures -- each folder can be utilized as a stand-alone analyses without accessing other data structures in this repo. To generate a figure, please run the associated Matlab script in each directory. Directories are divided by figure number and subsection. 


## LPS\_phenotyping<a id="lps_phenotyping"></a>

Raw gel images produced during O-antigen phenotyping are stored here (see DOI: [10.3791/3916](https://doi.org/10.3791/3916) for O-antigen extraction methods), as well as gel images from Lieberman and Michel’s 2011 work reanalyzed in this paper. 


## Supplemental\_figures<a id="supplemental_figures"></a>

Similar to “main\_figures”, the code needed to produce each supplementary figure is included here. Each folder can be utilized as a stand-alone analysis. In the event a figure requires sensitive health information, a note is included to explain that data’s absence.  


## Supplemental\_tables<a id="supplemental_tables"></a>

Scripts used to produce supplementary tables S3, S4, S5, S8, and S9 are included here, as well as a copy of the supplemental tables published in Poret et al. 
