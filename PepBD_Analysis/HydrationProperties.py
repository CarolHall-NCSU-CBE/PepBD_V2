#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 23 14:36:59 2021

@author: Michael
"""

import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy
import os
from pandas import read_csv
#from itertools import product

#dictionary relating 3letter amino acid names to 1 letter amino acid names
res_dict = {'ARG': 'R', "HIE": "H", "LYS" : "K", "ASP" : "D", "GLU" : "E",
            "SER" : "S", "THR": "T", "ASN" : "N", "GLN" : "Q", "CYS" : "C", 
            "GLY" : "G", "PRO" : "P", "ALA" : "A", "ILE" : "I", "LEU" : "L",
            "MET":"M", "PHE" : "F", "TRP": "W", "TYR": "Y", "VAL": "V"}

#KEY: 1=hydrophobic, 2=hydrophilic, 3=positive, 4=negative, 5=glycine, 6=other
types_long = {'ARG': 3, "HIE": 2, "LYS" : 3, "ASP" : 4, "GLU" : 4,
         "SER" : 2, "THR": 2, "ASN" : 2, "GLN" : 2, "CYS" : 6, 
         "GLY" : 5, "PRO" : 6, "ALA" : 6, "ILE" : 1, "LEU" : 1,
         "MET": 1, "PHE" : 1, "TRP": 1, "TYR": 1, "VAL": 1}

types_short = {'R': 3, "H": 2, "K" : 3, "D" : 4, "E" : 4,
         "S" : 2, "T": 2, "N" : 2, "Q" : 2, "C" : 6, 
         "G" : 5, "P" : 6, "A" : 6, "I" : 1, "L" : 1,
         "M": 1, "F" : 1, "W": 1, "Y": 1, "V": 1}

#masses of amino acids (in Daltons)
res_masses = {'R' : 174, "H" : 155, "K" : 146, "D" : 133, "E" : 147,
              "S" : 105, "T" : 119,"N" : 132, "Q" : 146, "C" : 121, 
              "G" : 75, "P" : 115, "A" : 89, "I" : 131, "L" : 131,
              "M" : 149, "F" : 165, "W" : 204,  "Y": 181, "V": 117}

### inputs used to run function if this script is called directly ###
plot_suffix='PepBD'
paths = ['/Users/Michael_1/']
file = 'TopPeps.csv'
hydration_single_run_flag = False #Change to true if you want the hydration profile for each run
hydration_net_average_flag = True #If true, then average hydration properties over all runs and plot
hydration_binary_flag = True #if true, then split amino acids into two categories: hydrophobic and hydrophilic

#KEY: 1=hydrophobic, 2=hydrophilic, 3=positive, 4=negative, 5=glycine, 6=other
type_labels_all=['Phobic', 'Philic', 'Pos', 'Neg', 'Gly', 'Other']
type_labels_binary = ['Polar', 'Non-Polar']

#plot parameters
#plt.style.use('/Users/Michael_1/Python/presentation.mplstyle') #this is my custom format file
colors = plt.get_cmap('Set1').colors
alpha=0.9

def read_pep_file(fullfile, file_type):
    try:
        if file_type == 'txt':
            with open(fullfile, 'r') as f:
                peps = f.readlines()
        elif file_type == 'csv':
            peps = read_csv(fullfile, header=0)['Seq']
        return peps
    except FileNotFoundError:
        return []

def calc_patch_size(pep):
    max_lengths = [0,0]
    cur_type = pep[0]
    length=1
    for AA in pep[1:]:
        if AA == cur_type:
            length += 1
        else:
            max_lengths[cur_type] = max( max_lengths[cur_type], length)
            cur_type = AA
            length = 1
    max_lengths[cur_type] = max( max_lengths[cur_type], length)
    return max_lengths

def HydrationProps_TopPeps(paths,
                           file,
                           hydration_single_run_flag=False,  
                           hydration_net_average_flag=False, 
                           verbose=False,
                           plot_suffix='',
                           plot_flag=True):

    #prepare storage for counting amino acid types
    pep_length = 0
    hist_density = True
    
    file_type = file.split('.')[-1]
    
    patch_phil = []
    patch_phob = []
    file_count = 0
    for path in paths:
        print(f'Analyzing {path}')

        peps = read_pep_file(os.path.join(path, file), file_type)
        if len(peps) == 0:
            print(f'Could not find information for {dir}! Skipping...')
            continue
        
        if pep_length == 0:
            pep_length = len(peps[0])
            net_counts_all = np.array([[0.0 for _ in range(6)] for _ in range(pep_length)])
            net_counts_binary = np.array([[0.0 for _ in range(2)] for _ in range(pep_length)])
        
        counts_all = np.array([[0.0 for _ in range(6)] for _ in range(pep_length)])
        counts_binary = np.array([[0.0 for _ in range(2)] for _ in range(pep_length)])
        for pep in peps:
            pep_all = [types_short[AA] for AA in pep]
            pep_binary = [0 if AA_type==1 or AA_type==6 else 1 for AA_type in pep_all]
            
            for i in range(pep_length):
                counts_all[i,pep_all[i]-1] +=1
                counts_binary[i,pep_binary[i]] +=1
                    
            #find largest patch size in each peptide (both hydrophobic and hydrophilic)
            patches = calc_patch_size(pep_binary)
            patch_phob.append(patches[0])
            patch_phil.append(patches[1])

        #calculate averages of amino acid counts
        net_counts_all = net_counts_all + counts_all/sum(counts_all[0])
        net_counts_binary = net_counts_binary + counts_binary/sum(counts_binary[0])
        
        if verbose:
            for i in range(len(type_labels_all)):
                print(f'Average Number of {type_labels_all[i]} amino acids: {sum(counts_all[:,i])}')
            print('\n')
            for i in range(len(type_labels_binary)):
                print(f'Average Number of {type_labels_binary[i]} amino acids: {sum(counts_binary[:,i])}')
            
        #plot results of single design run, if applicable
        if hydration_single_run_flag == True:
            pass
            
        file_count += 1

    #plot probability of finding each amino acid at each position in the peptide
    labels = [f'{i}' for i in range(1,pep_length+1)]
    net_counts_all = net_counts_all/file_count 
    net_counts_binary = net_counts_binary/file_count 
    if plot_flag == True:
        #first plot the results for all amino acid categories
        fig,ax=plt.subplots()
        ax.bar(labels, net_counts_all[:,0], label=type_labels_all[0], color=colors[0], alpha=alpha)
        boost = net_counts_all[:,0]
        for i in range(1,6):
            ax.bar(labels, net_counts_all[:,i], label=type_labels_all[i], bottom=boost, color=colors[i], alpha=alpha)
            boost += net_counts_all[:,i]
        ax.legend(bbox_to_anchor=(0.5,-0.4), bbox_transform=ax.transAxes, loc='center', ncol=2, fontsize=20)
        ax.set_xlabel('Peptide Residue')
        ax.set_ylabel('Frequency')
        ax.set_title(f'{plot_suffix}')
        
        #custom title
        #ax.set_title(f"CamSol, SF={d.split('_')[-1]}") #custom title
        
        #plot binary results
        fig,ax=plt.subplots()
        ax.bar(labels, net_counts_binary[:,0], label=type_labels_binary[0], color=colors[0], alpha=alpha)
        boost = net_counts_binary[:,0]
        for i in range(1,2):
            ax.bar(labels, net_counts_binary[:,i], label=type_labels_binary[i], bottom=boost, color=colors[i], alpha=alpha)
            boost += net_counts_binary[:,i]
        ax.legend(bbox_to_anchor=(0.5,-0.25), bbox_transform=ax.transAxes, loc='center', ncol=2, fontsize=20)
        ax.set_xlabel('Peptide Residue')
        ax.set_ylabel('Frequency')
        ax.set_title(f'{plot_suffix}')
               
    #make plot of patch size distributions of binary amino acid types
    fig,ax = plt.subplots()
    bins = np.arange(0,pep_length+1,1)    
    ax.hist(patch_phob, bins=bins, alpha=0.3, label='Hydrophobic', color=colors[0],  edgecolor='black', density=hist_density)
    ax.hist(patch_phil, bins=bins, alpha=0.3, label='Hydrophilic', color=colors[1],  edgecolor='black', density=hist_density)
    ax.legend(bbox_to_anchor=(0.5,-0.25), loc='center', bbox_transform=ax.transAxes, ncols=2, fontsize=20)
    ax.set_xlabel('Largest Patch Size')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{plot_suffix}')
    
if __name__ == '__main__': 
    HydrationProps_TopPeps(paths,
                           file,
                           hydration_single_run_flag, 
                           hydration_net_average_flag, 
                           plot_suffix=plot_suffix,
                           plot_flag=True)
    
    
    
    