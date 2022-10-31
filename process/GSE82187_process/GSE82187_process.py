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


data = 'GSE82187'

inpath = "./RAW/"
outpath = "../../%s/"%data
try:
    os.makedirs(outpath)
except:
    print()

#list of all the files

#meta data
meta = pd.read_csv(outpath + data + '_meta.csv')

#extract the sample name and sample title from the auxilary file
#Not sure why but meta data has duplicates with different  run name. GEO accession is the same with same cell type
#There are 2 types of protocol. I only sampled the ones with Mic-scRNA-Seq
aux = pd.read_csv(data + '_sample.csv')
raw_data = pd.read_csv(inpath + 'GSE82187_cast_all_forGEO.csv')



Exp_Name_list = []
Cell_type_list = []
for idx_p, protocol in enumerate(raw_data['protocol']):
    if protocol == 'Mic-scRNA-Seq':
        Exp_Name_list.append( raw_data['cell.name'][idx_p])
        Cell_type_list.append(raw_data['type'][idx_p])

Sample_Name_list = []
for exp_name in Exp_Name_list:
    found = False
    for idx_s, sample in enumerate(aux['Accession']):
        if exp_name in aux['Title'][idx_s]:
            found = True
            break
    if found == False:
        print('Exp name in raw data and aux file does not match')
    Sample_Name_list.append(sample)

no_dup = 0
for idx_s, sample in enumerate(Sample_Name_list):
    Exp_Name_list[idx_s] = Exp_Name_list[idx_s].split()[0]
    found = False
    for idx, geo in enumerate(meta['GEO_Accession (exp)']):
        if sample == geo:
            if found == False:
                found = True
                index = idx
            else:
                no_dup += 1
                if meta['major_cell_type'][idx] != meta['major_cell_type'][index]:
                    print('cell type of duplicate does not match!')
                    break
    if found == False:
        print('Sample name not found in meta data')



unique = list(set(Cell_type_list))
unique.sort()

N = len(Sample_Name_list)
print('Number of cells:', N)


#Do a quick check to make sure that raw data and labels match

cell_name_check = list(raw_data['cell.name'])
cell_type_check = list(raw_data['type'])
#make sure that the cell names match
if cell_name_check == Exp_Name_list:
    print('Cell name of raw data and meta matches!!')
if cell_type_check == Cell_type_list:
    print('Cell type of raw data and meta matches!!')


#construct the gene list
gene_list = list(raw_data.keys())[5:]   #the first row contains the gene info
for idx_gene, gene in enumerate(gene_list):
    gene_list[idx_gene] = gene.upper()

MATRIX = raw_data.values[:705, 5:].astype(float)
MATRIX = MATRIX.T

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
        
        

  
