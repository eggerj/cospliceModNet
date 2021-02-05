#!/usr/bin/env python3

'''
This program processes sample data from MAJIQ pipeline in order to create a sample by LSV PSI matrix and 
 LSV annotations for co-splicing network inference.
'''

# Import python libraries
import argparse
import operator
import numpy as np
import sys
sys.dont_write_bytecode = True

# Import custom libraries
#sys.path.insert(1, '/home/exacloud/lustre1/zheng_lab/users/eggerj/Dissertation/project_material/splicing_analysis/py_scripts')
import get_annotations as ga
import remove_redundant_LSVs as rr
import merge_lsv_clusters as mergeSVRs

# Load each sample file and collect LSV PSI values
def load_sample_data(fn):

    # Load in TSV file names
    with open(fn, 'r') as fh:
        sample_files = [l.strip() for l in fh.readlines()]

    # LSV dictionary of dictionaries
    #  {LSV_ID: {sample:line from tsv}}
    lsv_dict = {}

    # LSV dictionary of single record only for annotations
    lsv_recs = {}

    # Store sample names for later as well
    sample_names = []

    # Parse sample TSV files and collect LSV records
    for sf in sample_files:
        # Make generic
        sample = sf.split('/')[-1].split('.tsv')[0]
        sample_names.append(sample)
        # Load LSV records from current sample file
        with open(sf, 'r') as tsv_fh:
            lines = [l.strip().split('\t') for l in tsv_fh.readlines()]
        for rec in lines[1:]: # skip header
            lsv = rec[2] # get LSV ID
            # Add to LSV sample dictionary if it doesn't exist yet
            if lsv not in lsv_dict:
                lsv_dict[lsv] = {}
                lsv_recs[lsv] = rec
            lsv_dict[lsv][sample] = rec

    print('    ...total LSVs from sample set:', len(lsv_dict.keys()), '\n')

    return lsv_dict, sample_names, lsv_recs

# For unsupervised analysis (EDA and network inference), we need a single
#   value for each feature. Select junction having largest variance across
#    sample set.
def select_junctions(lsv_dict, lsv_recs):

    # Store selected PSI value, the junction used, and count of junctions in each LSV
    psi_dict = {}
    junction_indices = {}
    num_junctions = {}

    # For debugging
    #samples_index = [sample for samples in lsv_dict.values() for sample in samples.keys()][0]

    # Get numpy arrays of sample PSI values for every junction set of every LSV
    for lsv,samples in lsv_dict.items():
        # Create dictionary to store sample PSIs for each junction of LSV
        psi_lists = {j:[] for j in range(0,len(lsv_recs[lsv][3].split(';')))}
        for sample,rec in samples.items():
            junctionPSIs = rec[3].split(';')
            for j in range(len(junctionPSIs)):
                psi_lists[j] = psi_lists[j] + [float(junctionPSIs[j]),]
        # Determine which junction has greatest variance across samples
        psi_arrays = []
        for j in sorted(psi_lists.keys()):
            psi_arrays.append(np.array(psi_lists[j]))
        var_array = np.array([np.var(a) for a in psi_arrays])
        junc = np.argmax(var_array) # Index for which junction to use
        # Put LSV and sample PSIs in dictionary (dictionary of dictionaries --> {LSV:{sample:PSI}}
        # Store junction index with lsv as LSV_ID
        junction_indices[lsv] = str(junc + 1)
        num_junctions[lsv] = len(var_array)
        psi_dict[lsv] = {}
        # For each sample use PSI value from junction that had greatest variance
        for sample,rec in samples.items():
            junctionPSIs = rec[3].split(';')
            psi_dict[lsv][sample] = float(junctionPSIs[junc])
            
    return psi_dict, junction_indices, num_junctions

# Create sample by LSV PSI matrix and LSV annotations
def output_psi_matrix(psi_dict, junc_indices, redundant_lsvs, lsv_annotations, out_dir, group_name, num_juncs, lsv_dict, sample_names, lsv_regions, non_redundant):

    # Select samples and LSVs
    lsvs = sorted([lsv for lsv in psi_dict.keys()])
    samples = sorted(sample_names)
    sIndex = samples[0]

    fn = out_dir + '/psi_matrix.' + group_name + '.csv'

    # Create CSV file for holding PSI expression matrix of remaining LSVs
    with open(fn, 'w') as oh:
        oh.write(','+','.join(lsvs)+'\n')
        for sample in samples:
            sample_values = []
            for lsv in lsvs:
                if sample not in psi_dict[lsv]:
                    sample_values.append('NA')
                else:
                    sample_values.append(str(psi_dict[lsv][sample]))
            oh.write(sample + ',' + ','.join(sample_values) + '\n')

    # Create LSV annotations file for info and SVR formulation
    fn = out_dir + '/lsv_data_dictionary.' + group_name + '.csv'
    with open(fn, 'w') as oh:
        oh.write(',EVENT_TYPE,GENE_ID,GENE_NAME,LSV_TYPE,NUM_JUNCTIONS,JUNCTION,SVR,IS_REDUNDANT,REDUNDANT_PAIR\n')
        for lsv in lsvs:
            geneID = lsv_recs[lsv][1]
            geneName = lsv_recs[lsv][0]
            lsvType = lsv_recs[lsv][5]
            if lsv in non_redundant:
                nrLSV = 'keep'
            else:
                nrLSV = 'remove'
            oh.write(lsv + ',' + lsv_annotations[lsv] + ',' + geneID + ',' + geneName + ',' + lsvType + ',' + \
                     str(num_juncs[lsv]) + ',' + str(junc_indices[lsv]) + ',' + lsv_regions[lsv] + ',' + \
                     nrLSV + ',' + redundant_lsvs[lsv] + '\n') 


if __name__ == "__main__":

    # Program arguments
    parser = argparse.ArgumentParser(description='Create matrix for sample PSI values of LSVs and create file indicating SVR clusters (both CSV files).')
    requiredName = parser.add_argument_group('required arguments')
    requiredName.add_argument('--tsv_files', type=str, required = True, help='Text file containing paths and file names (one per line) of TSV files from voila tsv output.')
    requiredName.add_argument('--out_dir', type=str, required = True, help='Path to directory for writing output files.')
    requiredName.add_argument('--group_name', type=str, required = True, help='Identifier for dataset (will be file prefix).')
    requiredName.add_argument('--gtf', type=str, required = True, help='GTF file to get additional exon coordinate information.')
    args = parser.parse_args()

    # Get arguments and assign to variables
    fn = args.tsv_files
    out_dir = args.out_dir
    group_name = args.group_name
    gtf = args.gtf

    # Parse file of sample TSV file names and load LSVs 
    print('\n  Loading sample LSV data (' + group_name + ')...\n')
    lsv_dict, sample_names, lsv_recs = load_sample_data(fn)

    # Find LSVs that might be redundant (mirrors) and indicate in LSV dictionary output
    print ('  Finding potentially redundant LSVs...\n')
    sources, targets = rr.get_lsv_pairs(lsv_recs.values())
    redundant_pairs, mutually_exclusive = rr.find_redundant_pairs(sources, targets)
    lsv_dict, redundant_lsvs = rr.store_redundant_pairs(lsv_dict, redundant_pairs)
    non_redundant, x_hold = rr.filter_LSV_pairs(lsv_dict, redundant_pairs, mutually_exclusive) # Dont use x_hold
    print('    ...done.\n')

    # Annotate the event types of LSV set
    print('  Annotating LSVs...\n')
    # Load in exon coordinates from TSV files and GTF file
    gene_dict = ga.load_exon_coordinates(lsv_recs, gtf)    
    # Annotate LSVs
    lsv_annotations = ga.annotations(lsv_recs, gene_dict, mutually_exclusive)
    print('    ...done.\n')

    # Get overlapping LSVs and annotate as Splice Variant Regions (SVRs)
    print('  Assigning LSVs to SVRs...\n')
    lsv_regions = mergeSVRs.get_lsv_clusters(mergeSVRs.get_lsv_coords(lsv_recs))
    print('    ...done.\n')

    # Create expression matrix of LSVs (Rows = Samples, Columns = LSVs) and create table indicating 
    #  which LSVs should be merged into SVRs during processing (LSV dictionary).
    #  LSV dictionary also includes LSV features including gene ID/name, event type, and redundant pairs. 
    print('  Creating output files...\n')
    # Select which junction will be used for representing LSV in network (greatest variance across samples)
    psi_dict, junction_indices, num_juncs = select_junctions(lsv_dict, lsv_recs)
    # Create output
    output_psi_matrix(psi_dict, junction_indices, redundant_lsvs, lsv_annotations, out_dir, group_name, num_juncs, lsv_dict, sample_names, lsv_regions, non_redundant)
    print('    ...done.\n')


