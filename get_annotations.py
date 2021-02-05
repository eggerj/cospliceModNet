#!/usr/bin/env python3

'''
This library is used to create detailed event type annotations for LSVs formulated from MAJIQ splicing graphs.
'''

import argparse
import operator
import numpy as np
from collections import Counter


# Create merged set of exon coordinates
def merge_exon_coordinates(intervals):

    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)

    return merged


# This function gets all exon coordinates from both a GTF file and the
#  LSV splicing graph file from MAJIQ pipeline
def load_exon_coordinates(lsv_recs, gtf):

    # Create genes dictionary to store exon coordinates
    genes = {}

    # Add exon coordinates from LSV records
    for rec in lsv_recs.values():
        gene = rec[1]
        if gene not in genes:
            genes[gene] = []
        exons = rec[15].split(';')
        for exon in exons:
            if 'nan' not in exon.split('-'):
                exon_coords = tuple(sorted([int(exon.split('-')[0]), int(exon.split('-')[1])]))
                if exon_coords not in genes[gene]:
                    genes[gene] = genes[gene] + [exon_coords,]     

    # Now add additional exon coordinates from GTF file
    with open(gtf, 'r') as fh:
        lines = [l.strip().split('\t') for l in fh.readlines()]

    for rec in lines:
        gene = rec[8].split(';')[0].split('"')[1]
        if gene in genes:
            # get exon coordinates (3 and 4)
            exon_coords = tuple(sorted([int(rec[3]), int(rec[4])]))
            if exon_coords not in genes[gene] and rec[1] != 'retained_intron':
                genes[gene] = genes[gene] + [exon_coords,]

    # Sort exon coordinates for later
    for gene, exons in genes.items():
        genes[gene] = merge_exon_coordinates(sorted(genes[gene]))

    return genes

# Get number of exons covered by LSV junction
def get_exon_count(exons, junctions):

    max_junction = sorted(junctions, key = lambda sub: abs(sub[1] - sub[0]), reverse = True)[0]

    n = 0
    for exon in exons:
        if exon[0] >= max_junction[0] and exon[1] <= max_junction[1]:
            n = n + 1

    return n

def get_alt_splice_sites(junctions, rec_exons, gene_exons):

    # Get junction pairs resembling alt splice site usage
    # Junctions that share position on one end and use a different position in same exon on other end
    alt_pairs = []
    for j in range(0,len(junctions)-1):
        j1,j2 = junctions[j], junctions[j+1]
        if j1[0] != j2[0] and j1[1] == j2[1]:
            pair = sorted([j1[0], j2[0]])
        elif j1[0] == j2[0] and j1[1] != j2[1]:
            pair = sorted([j1[1], j2[1]])
        else:
            pair = None
        # See if junction pair shares exon at different sites
        if pair != None:
            alt_found = False
            for exon in rec_exons:
                # If pair shares exon it is alt splice site
                if pair[0] >= exon[0] and pair[1] <= exon[1]:
                    alt_found = True
                    alt_pairs.append((j1,j2))

    return alt_pairs

def get_exon_skip(junctions, rec_exons, gene_exons, alt_pairs, lsv):

    # Get junction pairs resembling exon skipping
    # Junctions that share one exon at one end and different exons at other end
    skipped_exons = 0
    exon_skip = False
    for j in range(0,len(junctions)-1):
        j1,j2 = junctions[j], junctions[j+1]
        if (j1,j2) not in alt_pairs:
            j1e1,j1e2, j2e1, j2e2 = None, None, None, None
            for exon in rec_exons:
                if j1[0] >= exon[0] and j1[0] <= exon[1]:
                    j1e1 = exon
                if j1[1] >= exon[0] and j1[1] <= exon[1]:
                    j1e2 = exon
                if j2[0] >= exon[0] and j2[0] <= exon[1]:
                    j2e1 = exon
                if j2[1] >= exon[0] and j2[1] <= exon[1]:
                    j2e2 = exon
            if (j1e1 == j2e1 and j1e2 != j2e2) or (j1e1 != j2e1 and j1e2 == j2e2):
                exon_skip = True
                n = get_exon_count(gene_exons, [j1,j2])
                if n > skipped_exons:
                    skipped_exons = n

    return exon_skip, skipped_exons

def annotations(lsvs, genes, mxe):

    exon_skips = []
    lsv_annotations = {}

    # Boolean annotations: 5'SS = 6, 3'SS = 7, ES = 8, IR = 16.split('|')[-1] (if == i)
    for lsv,rec in lsvs.items():
        alt5 = rec[6]
        alt3 = rec[7]
        es = rec[8]
        ir = rec[5].split('|')[-1]
        gene = rec[1]
        strand = rec[13]
        rec_exons = sorted([tuple(sorted([int(e.split('-')[0]), int(e.split('-')[1])])) for e in rec[15].split(';') if 'nan' not in e.split('-')])
        junctions = sorted([tuple(sorted([int(jr.split('-')[0]), int(jr.split('-')[1])])) for jr in rec[14].split(';')])
        lsv_type = lsv.split(':')[1]
        max_junction = sorted(junctions, key = lambda sub: abs(sub[1] - sub[0]), reverse = True)[0]
        # Mutually exclusive exons (found from remove_redundant_LSVs script)
        if lsv in mxe:
            lsv_annotations[lsv] = 'mutually_exclusive_exons'
        # Simple intron retentions
        elif es == 'False' and alt5 == 'False' and alt3 == 'False' and ir == 'i':
            lsv_annotations[lsv] = 'intron_retention'
        # Exclusively alternative splice sites
        elif es == 'False' and (alt5 == 'True' or alt3 == 'True') and ir != 'i':
            if len(junctions) == 2:
                lsv_annotations[lsv] = 'binary_altSS'
            else:
                lsv_annotations[lsv] = 'mult_altSS'     
        # Exclusively exon skipping events
        elif es == 'True' and alt5 == 'False' and alt3 == 'False' and ir != 'i':
            n = get_exon_count(genes[gene], junctions)
            if n > 0:
                exon_skips.append(n)
            if n <= 1:
                lsv_annotations[lsv] = 'exon_skip'
            else:
                lsv_annotations[lsv] = 'exon_cluster'
        # Intron retention with exon skipping
        elif es == 'True' and alt5 == 'False' and alt3 == 'False' and ir == 'i':
            n = get_exon_count(genes[gene], junctions)
            if n == 0:
                lsv_annotations[lsv] = 'intron_retention'
            else:
                exon_skips.append(n)
                lsv_annotations[lsv] = 'compound'
        # Alt splice site with exon skipping and/or IR
        elif (alt5 == 'True' or alt3 == 'True'): # and ir != 'i':
            # Check if alt splice site actually exists      
            alt_pairs = get_alt_splice_sites(junctions, rec_exons, genes[gene])
            num_alt_pairs = len(alt_pairs)
            # If not altSS, may still be exon skipping and/or IR
            if num_alt_pairs == 0: # and ir != 'i':
                if es == 'True':
                    n = get_exon_count(genes[gene], junctions)
                else:
                    n = 0
                if n == 0 and ir == 'i':
                    lsv_annotations[lsv] = 'intron_retention'
                elif n >= 1 and ir == 'i':
                    lsv_annotations[lsv] = 'compound'
                elif n == 1 and ir != 'i':
                    lsv_annotations[lsv] = 'exon_skip'
                elif n > 1 and ir != 'i':
                    lsv_annotations[lsv] = 'exon_cluster'
                else:
                    lsv_annotations[lsv] = 'unknown'
                if n >= 1:
                    exon_skips.append(n)
            # Just alt splice sites
            elif num_alt_pairs > 0 and es != 'True' and ir != 'i':
                if num_alt_pairs > 1:
                    lsv_annotations[lsv] = 'mult_altSS'
                else:
                    lsv_annotations[lsv] = 'binary_altSS'
            # Alt splice sites with exon skipping and/or IR
            elif num_alt_pairs > 0 and (es == 'True' or ir == 'i'):
                if es == 'True':
                    exon_skip, n = get_exon_skip(junctions, rec_exons, genes[gene], alt_pairs, lsv)
                    if (exon_skip):
                        lsv_annotations[lsv] = 'compound'
                        exon_skips.append(n)
                    elif ir == 'i':
                        lsv_annotations[lsv] = 'compound'
                    else:
                        if num_alt_pairs > 1:
                            lsv_annotations[lsv] = 'mult_altSS'
                        else:
                            lsv_annotations[lsv] = 'binary_altSS'
                elif ir == 'i':
                    lsv_annotations[lsv] = 'compound'
                else:
                    if num_alt_pairs > 1:
                        lsv_annotations[lsv] = 'mult_altSS'
                    else:
                        lsv_annotations[lsv] = 'binary_altSS'  
        
        else:
            lsv_annotations[lsv] = 'unknown'

    return lsv_annotations   
 

