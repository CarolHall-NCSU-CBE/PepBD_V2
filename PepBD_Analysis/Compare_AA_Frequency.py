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
pep_res = 10
pol = 'PE'
single_AA_plot = True
pair_AA_plot = False

## SPECIFY PATHS FOR COMPARISONS

# paths = [f'/Users/Michael_1/PeptideDesigns/PE/QuantumCircuit',
#          '/Users/Michael_1/Desktop/Temp/PepBD_Results/PE/MultipleStartStates']
# names = ['Quantum Circuit', 'PepBD']

# paths = [ '/Users/Michael_1/Desktop/Temp/PepBD_Results/PE/MultipleStartStates',
#          '/Users/Michael_1/Desktop/Temp/PepBD_Results/PS/MultipleStartStates',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PE_not_PS_01Dec2023']
# names = ['PepBD Designs: PE', 'PepBD Designs: PS', 'MCTS Designs: PE', 'MCTS Designs: PE vs. PS']

paths = [ '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianhong_Designs/PaperDesigns_13Oct2023/CamSol_2.0',
         '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PS/ML/Tianhong_Designs/PS_13Dec23',
         '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianhong_Designs/PE_not_PS_01Dec2023',
         '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PS/ML/Tianhong_Designs/PS_not_PE_11Dec23'
]
names = ['PE', 'PS', 'PE not PS', 'PS not PE']

paths = [ '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/QA+RL/QA_2Body_01Dec23',
         '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/QA+RL/QA+RL_GB+2Body_22Dec23',
         '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/PepBD/20Structs_FixedBB_PepBD']
names = ['QA', 'QA+RL', 'PepBD']

paths = [ '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PET/VQC_Raul',
         '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PepBD_Data/PepBD_Paper1/PET/SA+Hyd']
names = ['VQC', 'PepBD']

# paths = ['/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/QA+RL/QA+RL_GB+2Body_22Dec23',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/PepBD/20Structs_FixedBB_PepBD',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/QA+RL/RL_Match_PepBD_Score']
# names = ['QA+RL', 'PepBD', 'RL, Score Matching']

# paths = ['/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/CamSol_1.0',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/CamSol_2.0',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/CamSol_5.0',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/CamSol_10.0',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/CamSol_20.0']
# names = ['SF=1.0', 'SF=2.0', 'SF=5.0', 'SF=10.0', 'SF=20.0']

# paths = ['/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PepBD_Data/PepBD_Paper1/PE/SA+Hyd',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/Nominal',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/constrained']
# names = ['PepBD', 'MCTS', 'MCTS Constrained']

# paths = ['/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/ttH_0.01',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/ttH_0.1',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/ttH_0.3',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/ttH_0.5',
#          '/Users/Michael_1/Documents/Grad_School/HallGroup/MicroplasticsProject/PeptideDesign/PE/ML/Tianong_Designs/PaperDesigns_13Oct2023/ttH_1.0']
# names = ['SF=0.01', 'SF=0.1', 'SF=0.3', 'SF=0.5', 'SF=1.0']

paths = ['/Users/Michael_1/']
names = ['dummy']
num_names = len(names)

#plotting parameters
plt.rcdefaults()
#plt.style.use(f'{home}/Python/presentation.mplstyle') #this is my custom format file
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

if single_AA_plot:
    plot_single_AA_frequency(paths, pep_res)
if pair_AA_plot:
    plot_pair_AA_frequency(paths, pep_res)