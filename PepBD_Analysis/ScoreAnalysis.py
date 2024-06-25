#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 20:28:16 2022

@author: Michael
"""

from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import os
home = os.path.expanduser("~")

#file paths and names
file = 'energyprofile.txt'
rmsd_file = 'rmsd.txt'
out_file_1 = 'Score_Analysis_Results.csv'
out_file_2 = 'BestPeps_EnergyBreakdown.csv'

#analysis parameters
E_cut = 0 #if energy above this value, then don't include in analysis

#plot parameters

tick_step = 10
alpha=0.9
cur_cycle = -1
bin_width = 5
colors = plt.get_cmap('Pastel1').colors
linewidth=2
legend_fontsize=12
plt.style.use(os.path.join(home, 'Python/presentation.mplstyle')) #this is my custom format file
colors=plt.colormaps['Set1'].colors

#parameters used if calling function directly
paths = ['/Users/Michael_1/']
E_evolve_plot_flag = True #if True plot the evolution of PepBD score over course of design
rmsd_flag = True #if True, then plot evolution of RMSD as well as the score
E_evolve_multiplot_flag = False #if True, then plot evolution of all designs on a single plot
score_hist_flag = True #if True, then plot histogram of best scores from each design
score_hist_spline_fit_flag = True #if True, then fit spline to histogram
E_evolve_decomp_plot_flag = True #if True, then also plot the evolution of components of PepBD score during design
best_pep_flag = True #if True, then plot histogram of best scores over all designs

def ScoreAnalysis(paths,
                  E_evolve_plot_flag=True, 
                  E_evolve_multiplot_flag=True, 
                  score_hist_flag= False, 
                  score_hist_spline_fit_flag= False, 
                  rmsd_flag=True, 
                  E_evolve_decomp_plot_flag=False, 
                  best_pep_flag=True):
    
    if E_evolve_decomp_plot_flag == True and E_evolve_multiplot_flag == True:
        print('Both energy decomposition and multiple design runs should not be plotted simultaneously, the plots will be too difficult to read!')
        print('Turning off plotting of energy decomposition\n')
        E_evolve_decomp_plot_flag = False
    
    #prepare plots
    if score_hist_flag:
        fig2, ax2 = plt.subplots()
        ax2.set_xlabel('Total Binding Energy (kcal/mol)')
        ax2.set_ylabel('Probability')
        ax2.set_title('Histogram of PepBD Scores')
        
    if score_hist_spline_fit_flag:
        fig3, ax3 = plt.subplots()
        ax3.set_xlabel('Total Binding Energy (kcal/mol)')
        ax3.set_ylabel('Probability')
        ax3.set_title('PepBD Scores Distribution')
        
    if E_evolve_plot_flag == True and E_evolve_multiplot_flag == True:
        fig_evolve, ax_evolve = plt.subplots()
    c_index = 0
        
    #start of the actual analysis part
    best_scores = []
    all_scores = []
    out1 =  open(os.path.join(paths[0],out_file_1), 'w')
    out2 = open(os.path.join(paths[0],out_file_2), 'w')
    out1.write('Path,Score Mean,Score STD,Final RMSD\n')
    out2.write('Design,VDW,EGB,Peptide Energy,E_Bind\n')
    for path in paths:
            
        fullfile = os.path.join(path,file)
        count = 0
        try:
            #first figure out how many lines are in the file
            with open(fullfile, 'r') as f:
                all_data = f.read().splitlines()
                print(f'Processing file for {fullfile}')
        except FileNotFoundError:
            print(f'No file found for {os.path.join(path)}!')
            continue
             
        #now get the data for the file 
        count= len(all_data)
        step = np.empty(count)
        G_tot = np.empty(count)
        VDW = np.empty(count)
        ELE = np.empty(count)
        GB = np.empty(count)
        NP = np.empty(count)
        Pep_E = np.empty(count)
        G_trim = []
        index = 0
        best_E = E_cut
        for line in all_data[1:]:
            vals = line.split()
            if len(vals) < 8 :
                continue #this happens when energies are so large that there are not spaces between energies --> we don't care about this data
            vals[2:] = [float(v) for v in vals[2:]]
            step[index] = index
            G_tot[index] = vals[2]
            VDW[index] = vals[3]
            ELE[index] = vals[4]
            GB[index] = vals[5]
            NP[index] = vals[6]
            Pep_E[index] = vals[7]
            if G_tot[index] < E_cut:
                G_trim.append(G_tot[index])
            if G_tot[index] < best_E:
                best_index = index
                best_E = G_tot[index]
            index += 1
            
        #append minimum energy to array 
        best_scores.append(best_E)
        all_scores += list(G_trim)
                        
        #measure RMSD of peptide during design
        if rmsd_flag:
            try:
                data = open(os.path.join(path, rmsd_file)).readlines()
                rmsd = [float(line.split()[-1]) for line in data]
                final_rmsd = rmsd[-1]
            except FileNotFoundError:
                final_rmsd = 0.0
        else:
            final_rmsd = 0.0
                   
        av = np.average(G_trim)
        std = np.std(G_trim)
        print(f'{path}: mean = {av:2.2f}, std = {std:2.2f}, best = {best_E}, final RMSD: {final_rmsd}\n\n')
        out1.write(f'{path},{av:2.2f},{std:2.2f},{best_E:.2f},{final_rmsd:.2f}\n')

        #output energy decomposition of top peptide from run
        out2.write(f'{path},{VDW[best_index]},{ELE[best_index]},{NP[best_index]},{G_tot[best_index]}\n')
            
        #plot evolution of score for current design run
        if E_evolve_plot_flag == True:
            if E_evolve_multiplot_flag == False:
                fig_evolve, ax_evolve = plt.subplots()
            ax_evolve.plot(step[:index], G_tot[:index], label='Score', color=colors[c_index], alpha=alpha, linewidth=linewidth)
            if E_evolve_decomp_plot_flag == True:
                ax_evolve.plot(step[:index], VDW[:index], label='VDW', color=colors[c_index+1], alpha=alpha, linewidth=linewidth)
                ax_evolve.plot(step[:index], ELE[:index], label='ELE', color=colors[c_index+2], alpha=alpha, linewidth=linewidth)
                ax_evolve.plot(step[:index], GB[:index], label='GB', color=colors[c_index+3], alpha=alpha, linewidth=linewidth)
                ax_evolve.plot(step[:index], NP[:index], label='Peptide Energy', color=colors[c_index+4], alpha=alpha, linewidth=linewidth)
                ax_evolve.legend(bbox_to_anchor=(0.52, -.25), loc='center', ncol=5, fontsize=legend_fontsize)
                ax_evolve.set_ylim(min(min(ELE[:index]),min(VDW[:index])) - 10, max(GB[:index]) + 10)
            else:
                ax_evolve.set_ylim(best_E-10, E_cut)
            
            ax_evolve.set_title('PepBD Score Evolution')
            ax_evolve.set_xlabel('Step Number')
            if rmsd_flag == True:
                ax_r = ax_evolve.twinx()
                if E_evolve_decomp_plot_flag == True:
                    color='black'
                else:
                    color = colors[1]
                ax_r.plot(np.arange(1,len(rmsd)+1), rmsd, '-', color=color)
                ax_r.set_ylabel('RSMD (Å)', color=color)
                ax_evolve.set_ylabel('Energy (kcal/mol)', color=colors[0])
            else:
                ax_evolve.set_ylabel('Energy (kcal/mol)')

        #plot a histogram of the peptide scores for current design run
        if score_hist_flag:
            lim = [best_E, E_cut]
            nbins = int((E_cut - best_E)/bin_width)
            n, _, _ = ax2.hist(G_trim, bins=nbins, range=lim, density=True, alpha=0.5, label=path)
            
            #plot spline fit of the histogram heights (cubic spline)
            if score_hist_spline_fit_flag:
                nbins_interp = nbins*4
                interp_x_fit = np.linspace(lim[0], lim[1], num=nbins, endpoint=True)
                interp_x_plot = np.linspace(lim[0], interp_x_fit[-1], num=nbins_interp, endpoint=True)
                f = interp1d(interp_x_fit, n, kind='linear')
                ax3.plot(interp_x_plot, f(interp_x_plot), label=path)                                                 
        
    if score_hist_flag and len(paths)<10:
        ax2.legend(fontsize=legend_fontsize)
    if score_hist_spline_fit_flag and len(paths)<10:
        ax3.legend(fontsize=legend_fontsize)
        
    if best_pep_flag == True:
        fig,ax = plt.subplots()
        ax.set_title('Best PepBD Score Per Run')
        ax.set_xlabel('PepBD Score')
        ax.set_ylabel('Count')
        ax.hist(best_scores, color=colors[0], alpha=0.8, edgecolor='black')
        
    #make histogram of all scores from designs
    bin_min = np.floor(min(best_scores)/bin_width)*bin_width
    bin_max = np.ceil(max(best_scores)/bin_width)*bin_width
    bins = np.arange(bin_min, bin_max,bin_width, dtype=int)
    fig,ax = plt.subplots()
    ax.set_title('PepBD Score Distribution')
    ax.set_xlabel('PepBD Score')
    ax.set_ylabel('Count')
    ax.hist(all_scores, bins=bins, color=colors[0], alpha=0.8, edgecolor='black')
    
    out1.close()
    out2.close()

if __name__ == "__main__":
    ScoreAnalysis(paths,
                  E_evolve_plot_flag,
                  E_evolve_multiplot_flag,
                  score_hist_flag,
                  score_hist_spline_fit_flag,
                  rmsd_flag,
                  E_evolve_decomp_plot_flag,
                  best_pep_flag)
        
        
        
