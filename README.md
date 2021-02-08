# cospliceModNet


This repository contains Python and R scripts for performing co-splicing module network inference using complex splicing variants quantified from the MAJIQ framework for complex alternative splicing. Simulated local splicing variant (LSV) percent spliced in (PSI) values have been created in the form of output files from MAJIQ to use as an example. 

## Collecting LSV Splicing Values from Sample TSV Files

The script `create_psi_matrix.py` collects LSV splicing quantities from sample MAJIQ files (TSV files) and creates a .csv file representing a sample by LSV matrix of PSI quantities. In addition, the script creates a second .csv file containing annotations for each LSV including the gene name, gene ID (Ensembl), splicing event type, number of junctions, and most noteably, an SVR assignment indicating which SVR the LSV belongs to for SVR formulation (performed later in R). The remaining .py files contain functions used by `create_psi_matrix.py` for sample LSV processing. 

Run create_psi_matrix.py using Python 3 (`numpy` required):
```bash
python3 create_psi_matrix.py --tsv_files example_data/sample_LSVs.txt --out_dir test_out --gtf example_data/exons.example_data.gtf --group_name example_data  
```
Arguments:
- `tsv_files`: Text file listing names (and paths) of sample TSV files (one per line) created by MAJIQ (`majiq tsv`).
- `out_dir`: Name of directory (and path) for output files to be written.
- `gtf`: Exon models in GTF format for annotating event types of LSVs and SVRs. GTF should be the same build and version of the GFF file used during the build step with MAJIQ and should only contain "exon" records (third column). Exon only GTF can be created using awk.
- `group_name`: A unique identifier for file naming (e.g. name of experiment).

Example data is stored in the directory `example_data`:
- `sample_LSVs.txt`: Line by line listing of simulated sample LSV data. A total of 10,000 LSVs are present across the samples although some LSVs are missing in each sample (a common observation when setting read thresholds for PSI estimation with MAJIQ).
- `exons.example_data.gtf`: A GTF file containing exon records for all genes in which the simulated LSVs were derived. 

Output data created by `create_psi_matrix.py` using the example data has already been created and is stored in `test_out` for use with the R script `cosplicing_module_inference.R` (next section). 

## Formulating SVRs for De Novo Network Inference of Co-splicing Modules

The following R libraries are required for SVR formulation and network inference:
- `WGCNA`
- `stringr`  
- `ComplexHeatmap`
- `circlize`

The R script `cosplicing_module_infererence.R` loads the sample LSV quantities and LSV annotations in order to formulate SVRs and perform de novo network inference using the WGCNA framework.  


