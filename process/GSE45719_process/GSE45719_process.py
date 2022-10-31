# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 23:37:28 2022

@author: yutah
"""

import wget
import gzip, os
import numpy as np
import csv, os
import pandas as pd
from os import listdir
from os.path import isfile, join


data = 'GSE45719'



inpath = "./RAW/"
outpath = "../../%s/"%data
try:
    os.makedirs(outpath)
except:
    print()

try:
    os.makedirs(inpath)
except:
    print()



#list of all the files
onlyfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]   #the path of all the raw files

#meta data
meta = pd.read_csv(outpath + '%s_meta.csv'%(data))
source_name_list = meta['source_name']

GEO_list = list(meta['GEO_Accession (exp)'])
#construct the labels
Cell_type_list = []
Sample_Name_list = []
count = 0
for idx, source in enumerate(source_name_list):
    if '2-cell' in source:
        Cell_type_list.append('2-cell')
        Sample_Name_list.append(GEO_list[idx])
    elif '4-cell' in source:
        Cell_type_list.append('4-cell')
        Sample_Name_list.append(GEO_list[idx])
    elif '8-cell' in source:
        Cell_type_list.append('8-cell')
        Sample_Name_list.append(GEO_list[idx])
    elif '16-cell' in source:
        Cell_type_list.append('16-cell')
        Sample_Name_list.append(GEO_list[idx])
    elif 'Early blastocyst' in source:
        Cell_type_list.append('Early blastocyst')
        Sample_Name_list.append(GEO_list[idx])
    elif 'Mid blastocyst' in source:
        Cell_type_list.append('Mid blastocyst')
        Sample_Name_list.append(GEO_list[idx])
    elif 'Late blastocyst' in source:
        Cell_type_list.append('Late blastocyst')
        Sample_Name_list.append(GEO_list[idx])
    elif 'fibroblast' in source:
        Cell_type_list.append('fibroblast')
        Sample_Name_list.append(GEO_list[idx])
    else:
        print(source)
        count += 1

uq = list(set(list(Cell_type_list)))



N = len(Sample_Name_list)
print('Number of cells:', N)

#construct the gene list
gene_list = []
sample = Sample_Name_list[0]
for file in onlyfiles:
    if sample in file:
        found = True
        break
if not found:
    print(sample, 'not found')
raw_data_file = gzip.open(inpath + file, 'r')
raw_data_lines = raw_data_file.readlines()
raw_data_file.close()
gene_list = []
for idx_raw, raw_data in enumerate(raw_data_lines):
    if idx_raw > 0: #first line is header
        raw_data = raw_data.split()[0]  #first column is the gene name
        gene_list.append(raw_data.decode('utf-8').upper())

#Construct the data matrix
M = len(gene_list)
MATRIX = np.zeros([M, N])
for idx_cell in range(N):
    sample = Sample_Name_list[idx_cell]
    for file in onlyfiles:
        if sample in file:
            found = True
            break
    if not found:
        print(sample, 'not found')
    raw_data_file = gzip.open(inpath + file, 'r')
    raw_data_lines = raw_data_file.readlines()
    raw_data_file.close()
    for idx_raw, raw_data in enumerate(raw_data_lines):
        idx_gene = idx_raw - 1
        if idx_raw > 0: #first line is 
            raw_data = raw_data.split()
            gene = raw_data[0].decode('utf-8')
            if gene.upper() == gene_list[idx_gene]:
                MATRIX[idx_gene, idx_cell] = float(raw_data[2].decode('utf-8'))
            else:
                print('Gene Does not match')


#################
#################
#Write data for full data
#cell_type_unique = list(set(list(cell_type_full_list)))  #used if the cell_type makes sense
print('Writing the Full Data>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
cell_type_unique = ['2-cell', 'Late blastocyst', 'fibroblast', '8-cell', '16-cell', '4-cell', 'Early blastocyst', 'Mid blastocyst']

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

with open(outpath + data + '_full_labels.csv', 'w', newline='') as file:
    writer = csv.writer(file)
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
        
