# cospliceModNet


This repository contains Python and R scripts for performing co-splicing module network inference using complex splicing variants quantified from the MAJIQ framework for complex alternative splicing.

Simulated local splicing variant (LSV) percent spliced in (PSI) values have been created in the form of output files from MAJIQ to use as an example. 

The script create_psi_matrix.py collects LSV splicing quantities from sample MAJIQ files (TSV files) and creates a .csv file representing a sample by LSV PSI matrix. In addition, the script creates a second .csv file containing annotations for each LSV including the gene name, gene ID (Ensembl), splicing event type, number of junctions, and most noteably, an SVR assignment indicating which SVR the LSV belongs to for SVR formulation (performed later in R). The remaining .py files contain functions used by create_psi_matrix.py for sample LSV processing. 
