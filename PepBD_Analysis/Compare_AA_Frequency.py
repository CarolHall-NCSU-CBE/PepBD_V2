#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 14:49:35 2022

@author: Michael
"""

from matplotlib import pyplot as plt
import numpy as np
from pandas import read_csv
import os
home = os.path.expanduser("~")

#basic input files and parameters
file = 'AminoAcidFrequency.csv'
file_2 = 'AminoAcidPairFrequency.csv'
pep_res = 10 #specify number of amino acids per peptide
single_AA_plot = True #if True, then compare frequency of each amino acid
pair_AA_plot = False #if True, then compare frequency of amino acid dimers

# specify paths to amino acid frequency files, as well as names for the runs
paths = ['./']
names = ['PLACEHOLDER']
num_names = len(names)

#### STUFF BELOW HERE SHOULD NOT NORMALLY NEED TO BE CHANGED ####

#plotting parameters
bbox_loc = [0.5, -0.2]
colors = plt.get_cmap('Pastel2').colors
alpha=1
legend_fontsize=14
min_AA_space = 0.15
width = (1-min_AA_space)/num_names
shifter = num_names/2 - width

if num_names > 3:
    nrow=2
else:
    nrow=1
ncol = np.ceil(num_names/nrow)

### Analysis functions ###

def calc_enrichment(data, AAs, ref_i, data_i):
    #plot parameters
    min_AA_space = 0.15
    width = (1-min_AA_space)/num_names
    enrichment = np.log(data[data_i,:]/data[ref_i,:])
    fig, ax = plt.subplots(figsize=(10,5))
    x = np.arange(len(AAs))
    ax.set_xlabel('Amino Acid Type')
    ax.set_ylabel('Enrichment in PE Designs')
    ax.set_xticks(x, AAs)
    ax.hlines(y=0, xmin=x[0], xmax=x[-1], colors='black', linestyles='--')
    ax.bar(x, enrichment, width, alpha=alpha, color=colors[0], edgecolor='black')
    fig.tight_layout()

def plot_single_AA_frequency(paths, pep_res):
    #read and store all the data
    data = np.zeros((num_names,20))
    for index, path in enumerate(paths):
        df = read_csv(os.path.join(path, file))
        data[index,:] = df['Frequency']
    AAs = df['Amino_Acid']
               
    #prepare the plot
    fig, ax = plt.subplots(figsize=(10,5))
    x = np.arange(len(AAs))
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Frequency')
    ax.set_xticks(x, AAs)
    
    #plot the data
    for i in range(len(names)):
        if sum(data[i]) > 2:
            data[i] = data[i]/pep_res
        ax.bar(x + (i-shifter)*width, data[i], width, label=names[i], alpha=alpha, color=colors[i], edgecolor='black')
    
    ax.legend(bbox_to_anchor=bbox_loc, loc='upper center', ncol=ncol, fontsize=legend_fontsize)
    fig.tight_layout()
    
def plot_pair_AA_frequency(paths, cutoff_freq = 0.025):
    data = np.zeros((num_names,20**2))
    
    #read all data
    for i,name in enumerate(names):
        df = read_csv(os.path.join(paths[i], file_2))
        AAs = df['Amino_Acid_2mer']
        data[i,:] = df['Frequency'].to_numpy(dtype=float)
                
    #trim data set so only AA pairs with frequency above cutoff_freq are plotted
    AAs_trim= []
    data_trim = np.zeros((num_names,20**2))
    trim_index = 0
    for i in range(len(AAs)):
        if max(data[:,i] > cutoff_freq):
            AAs_trim.append(AAs[i])
            data_trim[:,trim_index] = data[:,i]
            trim_index += 1
    num_entries = len(AAs_trim)

    #do the plotting 
    fig, ax = plt.subplots(figsize=(10,5))
    x = np.arange(num_entries)
    ax.set_xlabel('Amino Acid 2mer')
    ax.set_ylabel('Frequency')
    ax.set_xticks(x, AAs_trim)
    for i in range(len(names)):
        ax.bar(x + (i-shifter)*width, data_trim[i,:num_entries], width, label=names[i], alpha=alpha, color=colors[i], edgecolor='black')
    
    ax.legend(bbox_to_anchor=bbox_loc, loc='upper center', ncol=ncol, fontsize=legend_fontsize)
    fig.tight_layout()
    
    
if __name__ == '__main__':
    if single_AA_plot:
        plot_single_AA_frequency(paths, pep_res)
    if pair_AA_plot:
        plot_pair_AA_frequency(paths, pep_res)