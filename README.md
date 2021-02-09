# cospliceModNet


This repository contains Python and R scripts for performing de novo network inference of co-splicing modules using complex alternative splicing variants. Complex splicing variants are annotated and quantified using the MAJIQ framework for complex alternative splicing. MAJIQ quantifies complex splicing in the form of local splicing variants (LSVs). Prior to network inference we formulate LSVs into splice variant regions (SVRs) representing variation in splicing levels across samples for a cluster of LSVs within a given gene. Simulated local splicing variant (LSV) percent spliced in (PSI) values have been created in the form of output files from MAJIQ to use as an example. 

## Creating Sample TSV Files of LSV PSI Quantities Using MAJIQ Framework

Using the MAJIQ framework we first create splicing graphs and annotate LSVs by running `majiq build` with all samples (.bam files) of our dataset. This will create, for each sample in our dataset, a file with the extension `.majiq` (e.g. `sample1.majiq`). These files are actually `.npz` files containing numpy arrays indicating the number of reads supporting each LSV junction in a given sample. 

In a typical differential splicing analysis we would first define our two sample groups and perform multiple hypothesis testing of relative splicing changes by running `majiq deltapsi`. For de novo network inference we need to estimate the relative abundance of each LSV across all samples of our dataset. Therefore, in order to obtain per-sample splicing values we run `majiq psi` on each sample `.majiq` file within our dataset. Running `majiq psi` on each sample will create a `.psi` file containing relative splicing values for each LSV meeting minimum thresholds in a given sample. These `.psi` files are to be used by the MAJIQ Voila functions. 

Under default settings, `majiq psi` will also create a `.tsv` alongside the `.psi` file. These are tab delimited files containing human readable splicing values for each LSV. For whatever reason, the `.tsv` files created here only contain limited information regarding each LSV (probably enough for most uses). For our needs, however, we require some more information regarding each LSV and we can obtain this information by generating more extensive `.tsv` files using MAJIQ's `voila tsv` function on each `.psi` file. This will create a `.tsv` with more columns than those created using `majiq psi` (we can set the `majiq psi` parameters to not generate the `.tsv` since we don't need them anyways).

For more details on how to work with the MAJIQ framework visit [majiq.biociphers.org](http://majiq.biociphers.org/)

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

Running `create_psi_matrix.py` will produce two output (`.csv`) files:
- `psi_matrix.example_data.csv`: PSI matrix where rows are samples and columns are LSVs
- `lsv_data_dictionary.example_data.csv`: LSV annotation dictionary where each row is an LSV and each column contains information regarding the LSV.


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

The R script `cosplicing_module_infererence.R` walks the user through loading the splicing data, formulating SVRs, and performing basic network inference and module detection using the WGCNA framework. It loads an additional file `cosplicing_helper_functions.R` which contains convenience functions for formulating SVRs and identifying co-splicing modules. Both output files from the previous step are required for SVR formulation.


## Credits & Acknowledgements
- [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/)
