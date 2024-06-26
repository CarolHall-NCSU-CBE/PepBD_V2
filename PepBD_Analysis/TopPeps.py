import numpy as np
import os

### PARAMETERS NEEDED TO RUN SCRIPT ###
length = 12 #number of amino acids in the peptide
top_sequence_numbers=50 #number of sequences per design run to extract
min_unique=3 #minimum number of amino acid differences needed in order to be considered a unique sequence

### PATHS AND FILE NAMES ###
paths = ['./'] #paths to where files are located, in list form
e_details_datafile = 'energydetails.txt' #name of file to analyze
outfile_all_name = 'TopPeps_ALL.csv' #name of file containing best sequences for all designs; arbitrarily created at first path location in "paths"

#### STUFF BELOW HERE SHOULD NOT NORMALLY NEED TO BE CHANGED ####

res_dict = {'ARG': 'R', "HIE": "H", "LYS" : "K", "ASP" : "D", "GLU" : "E",
            "SER" : "S", "THR": "T", "ASN" : "N", "GLN" : "Q", "CYS" : "C", 
            "GLY" : "G", "PRO" : "P", "ALA" : "A", "ILE" : "I", "LEU" : "L",
            "MET":"M", "PHE" : "F", "TRP": "W", "TYR": "Y", "VAL": "V"}

def check_bad_res(peptide, nonos):
    for res in peptide:
        if res in nonos:
            return 0
    return 1

def read_energydetails(e_details_file, decomp_flag=False):
    
    e_details_data = []
    pair = [0, 0, 0, 0]
    
    try:
        with open(e_details_file) as f:
            #determine features of data file
            info = f.readline().split()
            if info[2] == '1':
                subcycle_format_flag = 1
            else:
                subcycle_format_flag = 0
            f.readline()

            #collect data from the file
            f.seek(0)
            for line in f:
                info = line.split()
                if len(info) == 0:
                   continue
                elif info[0].isdigit():
                    if subcycle_format_flag == 0:
                        pair[0] = float(info[2])
                        pair[2] = int(info[0])
                        pair[3] = int(info[1])  
        
                    else:
                        pair[0] = float(info[3])
                        pair[2] = int(info[0])
                        pair[3] = int(info[1]) 
                            
                elif line[0] == 'N' or line[1]=='N':
                    info = line.split()
                    pair[1] = info
                    e_details_data.append(pair[:])
            
            #now sort the data to find top peptides with best binding free energy  
            temp = [elem[0] for elem in e_details_data]
            sort_inds = np.argsort(temp)
            
            return e_details_data, sort_inds
    
    except FileNotFoundError:
        return [], []

def is_pep_unique(new_seq, test_list, min_unique=3):
    for seq in test_list:
        if sum([a!=b for a,b, in zip(new_seq,seq)]) < min_unique:
            return False
    return True

def TopPeps(length, paths, top_sequence_numbers=50, min_unique=3, plot_flag=True):
    
    num_paths = len(paths)
    top_scores = np.zeros(num_paths) #stores top score from each design. Used to make a histogram of top scores at end
    for path_i, path in enumerate(paths):
        print (f'Processing {path}!')
        
        #read the energydetails file
        storage, sort_inds = read_energydetails(os.path.join(path, e_details_datafile))
        
        #check that the file could be found
        if len(storage) == 0:
            print (f'Did not find file in {path}! Skipping ...\n')
            continue
            
        #initialize variables
        cur_score = 0.0
        count = 0
        pep_list = []
        
        #for all directories, grab the top "target_num" peptides
        with open(os.path.join(path, 'TopPeps.csv'), 'w') as f_dir:
            f_dir.write('Score,Seq,Seq Long,Cycle,Step\n')
            for i in sort_inds:
                #grab next best scoring peptide
                [score, pep, cycle, step] = storage[i]
                
                #remove N and C from terminal amino acid residue names
                if len(pep[0]) == 4:
                    pep[0] = pep[0][1:]
                if len(pep[-1]) == 4:
                    pep[-1] = pep[-1][1:]
                
                if score!= cur_score: #simple check to make sure we don't count the same sequence twice
                    cur_score = score
                    
                    #convert peptide into one letter amino acid names
                    pep_long = ''
                    pep_short = ''
                    for res in pep:
                        pep_long += f'{res} '
                        pep_short += res_dict[res]
                        
                    #check if peptide is unique, given current list of peptides
                    if is_pep_unique(pep_short, pep_list):
                        pep_list.append(pep_short)

                        f_dir.write(f'{score:.2f},{pep_short},{pep_long},{cycle},{step}\n')
                    
                        count += 1
                        if count == top_sequence_numbers:
                            break
                        
    #write best_scoring peptide to net output file
    outfile_net = os.path.join(paths[0], outfile_all_name)
    with open(outfile_net,'w') as f:
        f.write('Path,Seq,Seq Long,Score\n')
        score, pep, cycle, step = storage[sort_inds[0]]
        pep_long = ''
        pep_short = ''
        for res in pep:
            pep_long += f'{res} '
            pep_short += res_dict[res]
        f.write(f'{path},{pep_short},{pep_long},{score}\n')
        top_scores[path_i] = score
            
    # #make histogram of best peptide scores (not currently used)
    # if len(top_scores)>1 and plot_flag == True:
    #     fig,ax = plt.subplots()
    #     ax.hist(top_scores, alpha=0.3,color='b', edgecolor='b')
    #     ax.set_xlabel('Top Peptide Score (kcal/mol)')
    #     ax.set_ylabel('Frequency')
    #     ax.set_title('Histogram of Best Peptide Scores')
    #     ax.set_title('Score Variability, 1 Conformation')

if __name__ == '__main__':

    TopPeps(length, paths,
            top_sequence_numbers=top_sequence_numbers, 
            min_unique=min_unique)