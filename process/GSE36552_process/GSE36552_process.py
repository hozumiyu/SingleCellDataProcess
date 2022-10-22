# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 10:50:23 2022

@author: yutah
"""

import gzip
import numpy as np
import csv
import pandas as pd
from os import listdir
from os.path import isfile, join
import gzip
import sys

data = 'GSE36552'


label = True
gene = True
X = True
inpath = "./%s_RAW/"%(data)
outpath = "../../%s/"%data 

meta = pd.read_csv(outpath + '%s_meta.csv'%data)
onlyfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]   #the path of all the raw files
Accession_list = meta['GEO_Accession (exp)']      #the sample name  

            
if label:
    #generate labels
    GEO_Accession = meta['GEO_Accession (exp)']
    cell_type_full_list = meta['Library Name']   #this line changes depending on the problem
    cell_type_unique = ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', 
                        'Morulae', 'Late blastocyst', 'hESC passage#0', 'hESC passage#10']
    #cell_type_unique = ['Zygote', '2-cell', '4-cell', '8-cell', '16cell', 'Late blastocyst', 'Oocyte', 'Morulae']
    #cell_type_unique = list(set(list(cell_type_full_list)))  #used if the cell_type makes sense
    cell_type_dict = {i: cell_type_unique[i] for i in range(len(cell_type_unique)) }  #dict to map index to cell type
    cell_type_dict_rev = {cell_type_unique[i]:i for i in range(len(cell_type_unique))}  #reverse dict from cell type to index
    file = open(outpath + data + '_labeldict.txt', 'w')   #write a dictionary
    file.writelines(str(cell_type_dict) + ' \n')
    file.writelines(str(cell_type_dict_rev) + ' \n')
    file.close()
    with open(outpath + data + '_labels.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sample Name', 'Cell type', 'Label'])
        for idx in range(len(cell_type_full_list)):
            cell_type = cell_type_full_list[idx]
            name = None
            for c in cell_type_unique:
                if c in cell_type:
                    name = c
                    
                    break
            #if name == 'Morulae':
            #    name = '16cell'
            #elif name == 'Oocyte':
            #    name = 'Zygote'
            if name == None:
                if 'ES p0#1' in cell_type:
                    name = 'hESC passage#0'
            if name == None:
                print('Error')
                break
            if name:
                writer.writerow([GEO_Accession[idx], name, cell_type_dict_rev[name]])
            #if name == None:
            #    if 'ES p0#1' in cell_type:
            #        name = 'hESC passage#0'
            #if name == None:
            #    print('Error')
            #    break
            
            


'''
this is used when the gene is different between cells
'''
if gene:  
    sample_name_list = pd.read_csv(outpath + data + '_labels.csv' )['Sample Name']
    #generate the genes
    full_gene = []
    for sample in sample_name_list:
        for file in onlyfiles: # get the file name associated with sample name
            if sample in file:
                break
        file = gzip.open(inpath + file, 'r')  # get the file
        raw_data_lines = file.readlines()
        file.close()
        for idx_line, line in enumerate(raw_data_lines):
            if idx_line > 0:        #skip the first line
                gene = line.split()[0].decode('utf-8') #first index is the gene
                if gene not in full_gene:  #only keep the genes that was not defined earlier
                    full_gene.append(gene)
       
    #full_gene.sort()
    with open(outpath + data + '_gene.csv', "w", newline = '') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['Index', 'Gene'])
        for idx, gene in enumerate(full_gene):
            writer.writerow([idx, gene])
    print('Number of gene:', len(full_gene))
    

if X:
    gene_list = list(pd.read_csv(outpath + data + '_gene.csv')['Gene']   )#gene list that will be used to compare the files
    sample_name_list = list(pd.read_csv(outpath + data + '_labels.csv' )['Sample Name'])
    Matrix = np.zeros([len(gene_list), len(sample_name_list)])
    for idx_col, sample in enumerate(sample_name_list):
        print(idx_col)
        found = False
        for file in onlyfiles: # get the file name associated with sample namea
            if sample in file:
                found = True
                break
        if found == False:
            print(sample, 'file does not exist')
        file = gzip.open(inpath + file, 'r')  # get the file
        raw_data_lines = file.readlines()
        file.close()
        for idx_line, line in enumerate(raw_data_lines):
            if idx_line > 0:
                line = line.split()
                gene = line[0].decode('utf-8')
                rpkm = eval(line[-1])
                for idx_row, g in enumerate(gene_list):  #get the row associated with the gene
                    found = False
                    if g == gene:
                        found = True
                        Matrix[idx_row, idx_col] = rpkm
                        break
                if found == False:
                    print('gene not found', print(gene))
    
    with open(outpath + data + '_data.csv', "w", newline = '') as csv_file: #output file for data
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['Row', 'Col', 'Val'])
        for idx_col in range(len(sample_name_list)):
            for idx_row in range(len(gene_list)):
                if Matrix[idx_row, idx_col] > 0:
                    writer.writerow([idx_row, idx_col, Matrix[idx_row, idx_col]])
    #Matrix[Matrix < 0.1] = 0
    with open(outpath + data + '_X.csv', "w", newline = '')  as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow([None] + sample_name_list)
        for idx_row in range(len(gene_list)):
            row = list(Matrix[idx_row, :])
            writer.writerow([gene_list[idx_row] ] + row)
    
