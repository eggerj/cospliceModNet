#!/usr/bin/env python3

'''
This program is for identifying potential redundant LSVs that could be removed prior to SVR formulation.
  Redundant LSVs can be optionally removed, but are usually kept for reference (can be removed later).
'''


import argparse
import operator
import numpy as np


def get_shared_LSVs(lsv_dict, samples, min_sample_ratio):

    shared_lsvs = {}

    for lsv,data in lsv_dict.items():
        lsv_ratio = float(len(data.keys())) / float(len(samples))
        if lsv_ratio < min_sample_ratio:
            foo = 1 #TODO: write missing LSVs to file
        else:
            shared_lsvs[lsv] = data

    print('Total LSVs:', len(lsv_dict.keys()))
    print('Total shared LSVs:', len(shared_lsvs.keys()))

    return shared_lsvs

def get_lsv_pairs(lsv_set):

    # Dictionary of dictionaries: {geneID: {lsv_id:rec}}
    source_LSVs = {}
    target_LSVs = {}

    # Collect LSVs that do not contain alt splice sites or intron retentions
    for rec in lsv_set:
        lsv_type = rec[5]
        alt5 = rec[6]
        alt3 = rec[7]        
        if lsv_type.split('|')[-1] != 'i' and alt5 != 'True' and alt3 != 'True':
            gene_id = rec[1]
            lsv_id = rec[2]
            if lsv_id.split(':')[1] == 's':
                if gene_id not in source_LSVs:
                    source_LSVs[gene_id] = {}
                source_LSVs[gene_id][lsv_id] = rec
            else: 
                if gene_id not in target_LSVs:
                    target_LSVs[gene_id] = {}
                target_LSVs[gene_id][lsv_id] = rec

    return source_LSVs, target_LSVs       

def search_rec(source_rec, target_rec):

    # Look for the following:
    #  1) Share a junction
    #  2) Share first and last exons if 2 junctions exist
    #  3) Share all exons if 3 or more junctions exist 

    sourceID = source_rec[2]
    targetID = target_rec[2]
    # Pull junctions from records
    source_juncs = source_rec[14].split(';')
    target_juncs = target_rec[14].split(';')
    # Pull exons from records
    source_exons = sorted([sorted(pair.split('-')) for pair in source_rec[15].split(';')])
    target_exons = sorted([sorted(pair.split('-')) for pair in target_rec[15].split(';')])

    junction_shared = False
    for sj in source_juncs:
        if sj in target_juncs:
            junction_shared = True
            break

    if (junction_shared) and len(source_juncs) == len(target_juncs):
        # If binary LSVs, only first and last exons have to match
        if len(source_juncs) == 2:
            if source_exons[0] == target_exons[0] and source_exons[-1] == target_exons[-1]:
                return True, False
        # If 3+ junctions, all exons have to match
        else:
            if source_exons == target_exons:
                return True, False
    # Find mutually exclusive exons?
    else:
        if len(source_exons) == 3 and len(target_exons) == 3 and len(source_juncs) == 2 and len(target_juncs) == 2:
            if source_exons[1] == target_exons[0] and source_exons[2] == target_exons[1]:
                return True, True

    return False, False

def find_redundant_pairs(sources, targets):

    redundant_pairs = []
    mutually_exclusive = {}

    for gene_id, source_LSVs in sources.items():
        if gene_id in targets:
            for source_lsv, source_rec in source_LSVs.items(): # source LSVs for current gene
                for target_lsv, target_rec in targets[gene_id].items(): # target LSVs for current gene
                    pair_found, mxe = search_rec(source_rec, target_rec)
                    if (pair_found):
                        redundant_pairs.append((source_lsv, target_lsv, mxe))
                    if (mxe):
                        mutually_exclusive[source_lsv] = 'mutually_exclusive'
                        mutually_exclusive[target_lsv] = 'mutually_exclusive'

    return redundant_pairs, mutually_exclusive

def get_variance(lsv, lsv_dict):

    junc_vars = [v[4].split(';') for v in lsv_dict[lsv].values()]
    return np.array([float(v) for var in junc_vars for v in var])

def filter_LSV_pairs(lsv_dict, pairs, mut_excl):

    redundant_lsvs = {}
    filtered_lsvs = {}
    removed_lsvs = {}
    hold_pairs = []

    # Try this: Keep LSV having smallest average variance of expected PSI values across samples (average of junctions, average of samples)
    for pair in pairs:
        source,target,mxe = pair[0],pair[1],pair[2]
        # Check that one LSV is part of a MXE set and if so save for later
        if (mxe == False) and (source in mut_excl or target in mut_excl):
            hold_pairs.append(pair)
        # If source or target already removed, keep it that way and keep the other LSV 
        elif target in removed_lsvs:
            filtered_lsvs[source] = lsv_dict[source]  
            if source in redundant_lsvs:
                redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
            else:
                redundant_lsvs[source] = target
        elif source in removed_lsvs:
            filtered_lsvs[target] = lsv_dict[target]
            if target in redundant_lsvs:
                redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
            else:
                redundant_lsvs[target] = source
        # Else, select LSV with lowest variance between the two across samples
        elif np.mean(get_variance(source, lsv_dict)) > np.mean(get_variance(target, lsv_dict)): # Keep target LSV
            filtered_lsvs[target] = lsv_dict[target]
            removed_lsvs[source] = source
            if target in redundant_lsvs:
                redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
            else:
                redundant_lsvs[target] = source
        else: 
            filtered_lsvs[source] = lsv_dict[source]
            removed_lsvs[target] = target
            if source in redundant_lsvs:
                redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
            else:
                redundant_lsvs[source] = target

    # Take care of redundant LSVs inside of mutually exclusive exon pair
    for pair in hold_pairs:
        source,target = pair[0],pair[1]
        # Keep LSV that was kept to represent MXE and remove other LSV
        if source in mut_excl and source in filtered_lsvs:
            removed_lsvs[target] = target
            redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
        elif source in mut_excl and source not in filtered_lsvs:
            filtered_lsvs[target] = lsv_dict[target]
            removed_lsvs[source] = source
            if target in redundant_lsvs:
                redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
            else:
                redundant_lsvs[target] = source 
        elif target in mut_excl and target in filtered_lsvs:
            removed_lsvs[source] = source
            redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
        else:
            filtered_lsvs[source] = lsv_dict[source]
            removed_lsvs[target] = target
            if source in redundant_lsvs:
                redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
            else:
                redundant_lsvs[source] = target

    # Store remaining LSVs redundant pairs as 'NaN' and store remaining LSVs without a redundant pairing 
    for lsv in lsv_dict.keys():
        if lsv not in redundant_lsvs and lsv not in removed_lsvs:
            redundant_lsvs[lsv] = 'NaN'
            filtered_lsvs[lsv] = lsv_dict[lsv] 

    return filtered_lsvs, redundant_lsvs

def store_redundant_pairs(lsv_dict, pairs):

    redundant_lsvs = {}

    for pair in pairs:
        source, target = pair[0], pair[1]
        if target in redundant_lsvs:
            redundant_lsvs[target] = redundant_lsvs[target] + '|' + source
        else:
            redundant_lsvs[target] = source
        if source in redundant_lsvs:
            redundant_lsvs[source] = redundant_lsvs[source] + '|' + target
        else:
            redundant_lsvs[source] = target 

    # Store remaining LSVs as having no redundant pairing
    for lsv in lsv_dict.keys():
        if lsv not in redundant_lsvs:
            redundant_lsvs[lsv] = 'NaN'    

    return lsv_dict, redundant_lsvs       


    
