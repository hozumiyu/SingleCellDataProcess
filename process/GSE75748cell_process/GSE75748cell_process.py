# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 23:37:28 2022

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


data = 'GSE75748cell'


inpath = "./RAW/"
outpath = "../../%s/"%data
try:
    os.makedirs(outpath)
except:
    print()

#list of all the files

#meta data
meta = pd.read_csv(outpath + data + '_meta.csv')
aux = pd.read_csv(data + '_sample.csv')
raw_data = pd.read_csv(inpath + 'GSE75748_sc_cell_type_ec.csv')


#header of the raw data cotains the experiment names
Exp_Name_list = list(raw_data.keys())[1:]


Sample_Name_list = []
Cell_type_list = [] 
for exp_name in Exp_Name_list:
    found = False
    for idx_s, sample in enumerate(aux['Accession']):
        title = aux['Title'][idx_s]
        title = title.split()[-1]
        title = title[1:-1]
        if exp_name == title:
            found = True
            break
    if found == False:
        print('Exp name in raw data and aux file does not match')
        break
    Sample_Name_list.append(sample)
    cell_type = exp_name.split('_')[0]
    Cell_type_list.append(cell_type)

#check to see if there is any sample name duplicates
if len(list(set(Sample_Name_list)))!= len(Exp_Name_list):
    print('Number of exp and sample does not match')

#construct the gene list
gene_list = list(raw_data['Unnamed: 0'])
for idx_gene, gene in enumerate(gene_list):
    gene_list[idx_gene] = gene.upper()


unique = list(set(Cell_type_list))
unique.sort()

N = len(Sample_Name_list)
M = len(gene_list)
print('Number of cells:', N)
print('Number of Genes:', M)

if len(raw_data) != M:
    print('Number of gene does not match raw data')



MATRIX = raw_data.values[:, 1:].astype(float)

if MATRIX.shape[0] == len(gene_list) and MATRIX.shape[1] == N:
    print('Dimension of the gene and cell matches!')
else:
    print('WARNING: DIMENSION OF THE MATRIX DOES NOT MATCH')

#################
#################
#Write data for full data
#cell_type_unique = list(set(list(cell_type_full_list)))  #used if the cell_type makes sense
print('Writing the Full Data>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
cell_type_unique = unique

#write the dictionary
cell_type_dict = {i: cell_type_unique[i] for i in range(len(cell_type_unique)) }  #dict to map index to cell type
cell_type_dict_rev = {cell_type_unique[i]:i for i in range(len(cell_type_unique))}  #reverse dict from cell type to index
file = open(outpath + data + '_full_labeldict.txt', 'w')   #write a dictionary
file.writelines(str(cell_type_dict) + ' \n')
file.writelines(str(cell_type_dict_rev) + ' \n')
file.close()
cell_count =  {cell_type_unique[i]:0 for i in range(len(cell_type_unique))}
Sample_Name_list = Sample_Name_list
Cell_label_list = []
for idx in range(len(Cell_type_list)):  #iterate over all the samples
    cell_type = Cell_type_list[idx]     #get the current cell type
    try:   #see if the cell type exist in the dictionary
        cell_label = cell_type_dict_rev[cell_type]
        cell_count[cell_type] += 1
        Cell_label_list.append(cell_label)
        found = True
    except:
        found = False    #if it is not found, print message
        print(Sample_Name_list[idx], cell_type, 'Wrong cell type')

print('Cell Count')
for k in cell_count.keys():
    print('%s:'%k, cell_count[k])

with open(outpath + data + '_full_labels.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Sample Name', 'Cell type', 'Label'])
    for idx in range(len(Cell_type_list)):
        writer.writerow([Sample_Name_list[idx], Cell_type_list[idx] , Cell_label_list[idx]])
        
                
with open(outpath + data + '_full_gene.csv', "w", newline = '') as csv_file:
    writer = csv.writer(csv_file, delimiter=',')
    writer.writerow(['Index', 'Gene'])
    for idx, gene in enumerate(gene_list):
        writer.writerow([idx, gene])
                


with open(outpath + data + '_full_data.csv', "w", newline = '') as csv_file: #output file for data
    writer = csv.writer(csv_file, delimiter=',')
    writer.writerow(['Row', 'Col', 'Val'])
    for idx_col in range(len(Sample_Name_list)):
        for idx_row in range(len(gene_list)):
            if MATRIX[idx_row, idx_col] > 0:
                writer.writerow([idx_row, idx_col, MATRIX[idx_row, idx_col]])
                    
with open(outpath + data + '_full_X.csv', "w", newline = '') as csv_file: #output file for data
    writer = csv.writer(csv_file, delimiter=',')
    writer.writerow( [None] + Sample_Name_list)
    for idx_row in range(len(gene_list)):
        row = list(MATRIX[idx_row, :])
        writer.writerow([gene_list[idx_row] ] + row)
        
        

  
