# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 10:50:23 2022

@author: yutah
"""

import gzip
import numpy as np
import csv, os
import pandas as pd
from os import listdir
from os.path import isfile, join
import gzip
import sys


data = 'GSE45719'

meta = False
label = True
gene = True

X = True
inpath = "./RAW/"
outpath = "../../%s/"%data
try:
    os.makedirs(outpath)
except:
    print()

onlyfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]   #the path of all the raw files

if label:
    #generate labels
    meta = pd.read_csv(outpath + data + '_meta.csv')  
    GEO_Accession = meta['GEO_Accession (exp)']
    cell_type_full_list = meta['source_name']   #this line changes depending on the problem
    # I removed the zygote, liver cells, fibroblast and 2-cell stage blastomere. Total of 9 classes
    cell_type_unique = [#'Zygote', 
                        'Early 2-cell stage blastomere (31-32h post-fertilization)',
                        'Mid 2-cell stage blastomere (34-40h post-fertilization)', 
                        'Late 2-cell stage blastomere (46-48h post-fertilization)',
                        '4-cell stage blastomere',
                        '8-cell stage blastomere',
                        '16-cell stage blastomere',
                        'Early blastocyst cell (86-88h post-fertilization)',
                        'Mid blastocyst cell (92-94h post-fertilization)',
                        'Late blastocyst cell (100-102h post-fertilization)']
    #cell_type_unique = list(set(list(cell_type_full_list)))  #used if the cell_type makes sense
    cell_type_dict = {i: cell_type_unique[i] for i in range(len(cell_type_unique)) }  #dict to map index to cell type
    cell_type_dict_rev = {cell_type_unique[i]:i for i in range(len(cell_type_unique))}  #reverse dict from cell type to index
    file = open(outpath + data + '_labeldict.txt', 'w')   #write a dictionary
    file.writelines(str(cell_type_dict) + ' \n')
    file.writelines(str(cell_type_dict_rev) + ' \n')
    file.close()
    cell_count =  {cell_type_unique[i]:0 for i in range(len(cell_type_unique))}
    Sample_Name_list = []
    Cell_type_list = []
    Cell_label_list = []
    for idx in range(len(cell_type_full_list)):  #iterate over all the samples
        cell_type = cell_type_full_list[idx]     #get the current cell type
        try:   #see if the cell type exist in the dictionary
            cell_label = cell_type_dict_rev[cell_type]
            cell_count[cell_type] += 1
            Cell_type_list.append(cell_type); Cell_label_list.append(cell_label); Sample_Name_list.append(GEO_Accession[idx])
            found = True
        except:
            found = False    #if it is not found, print message
            print(GEO_Accession[idx], cell_type, 'Wrong cell type')
    print('Number of Cells:', len(Cell_type_list))
    print(cell_count)
    
    with open(outpath + data + '_labels.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sample Name', 'Cell type', 'Label'])
        for idx in range(len(Cell_type_list)):
            writer.writerow([Sample_Name_list[idx], Cell_type_list[idx] , Cell_label_list[idx]])
        
                

    
    
'''
Use this if the gene is the same across samples
'''

if gene:
    #generate the genes
    labels_file = pd.read_csv(outpath + data + '_labels.csv')
    sample = labels_file['Sample Name'][0]
    found = False
    for file in onlyfiles:
        if sample in file:
            found = True
            break
    if found == False:
        print('File not found')
    
    raw_data_flie = gzip.open(inpath + file, 'r')
    raw_data_lines = raw_data_flie.readlines()
    raw_data_flie.close()
    gene_list = []
    for idx_raw, raw_data in enumerate(raw_data_lines):
        if idx_raw > 0: #first line is 
            raw_data = raw_data.split()[0]  #first column is the gene name
            gene_list.append(raw_data.decode('utf-8'))
    with open(outpath + data + '_gene.csv', "w", newline = '') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['Index', 'Gene'])
        for idx, gene in enumerate(gene_list):
            writer.writerow([idx, gene])
                
if X:
    labels_file = pd.read_csv(outpath + data + '_labels.csv')
    Sample_Name_list = list(labels_file['Sample Name'])
    gene_list = list(pd.read_csv(outpath + data + '_gene.csv')['Gene'])
    Matrix_rpkm = np.zeros([len(gene_list), len(Sample_Name_list)])
    Matrix_count = np.zeros([len(gene_list), len(Sample_Name_list)])
    for idx_col, sample in enumerate(Sample_Name_list):
        found = False
        for file in onlyfiles:
            if sample in file:
                found = True
                break
        if found == False:
            print(sample, 'file not found')
        #get the file
        raw_data_flie = gzip.open(inpath + file, 'r')
        raw_data_lines = raw_data_flie.readlines()
        raw_data_flie.close()
        for idx_raw, raw_data in enumerate(raw_data_lines):
            idx_row = idx_raw - 1
            if idx_raw > 0: #first line is 
                raw_data = raw_data.split()  #first column is the gene name
                gene = raw_data[0].decode('utf-8')
                if gene == gene_list[idx_row]: #check to make sure that the gene matches
                    rpkm = float(raw_data[2].decode('utf-8')) ; val = float(raw_data[3].decode('utf-8'))
                    if rpkm > 0:
                        Matrix_rpkm[idx_row, idx_col] = rpkm
                    if val > 0:
                        Matrix_count[idx_row, idx_col] = val
            
    with open(outpath + data + '_data.csv', "w", newline = '') as csv_file: #output file for data
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['Row', 'Col', 'Val'])
        for idx_col in range(len(Sample_Name_list)):
            for idx_row in range(len(gene_list)):
                if Matrix_rpkm[idx_row, idx_col] > 0:
                    writer.writerow([idx_row, idx_col, Matrix_rpkm[idx_row, idx_col]])
                    
    with open(outpath + data + '_X.csv', "w", newline = '') as csv_file: #output file for data
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow( [None] + Sample_Name_list)
        for idx_row in range(len(gene_list)):
            row = list(Matrix_rpkm[idx_row, :])
            writer.writerow([gene_list[idx_row] ] + row)
                
    with open(outpath + data + 'count_data.csv', "w", newline = '') as csv_file: #output file for data
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['Row', 'Col', 'Val'])
        for idx_col in range(len(Sample_Name_list)):
            for idx_row in range(len(gene_list)):
                if Matrix_count[idx_row, idx_col] > 0:
                    writer.writerow([idx_row, idx_col, Matrix_count[idx_row, idx_col]])
                    
    with open(outpath + data + 'count_X.csv', "w", newline = '') as csv_file: #output file for data
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow( [None] + Sample_Name_list)
        for idx_row in range(len(gene_list)):
            row = list(Matrix_count[idx_row, :])
            writer.writerow([gene_list[idx_row] ] + row)
    
