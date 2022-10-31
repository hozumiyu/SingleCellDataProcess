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


data = 'GSE63818'

label = True
gene = False

X = False
inpath = "./RAW/"
outpath = "../../%s/"%data
try:
    os.makedirs(outpath)
except:
    print()

#list of all the files
onlyfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]   #the path of all the raw files

#meta data
aux = pd.read_csv('GSE63818_sample.csv')
title_list = aux['Title']

#construct the labels
Cell_type_list = []
for title in title_list:
    title = title.split('_')
    title2 = '_'.join(title[:3])
    if 'F_PGC_8W' == title2:
        title[-1] = title[-1][:2]
        title2= '_'.join(title[:5])
        #print(title2)
        
    Cell_type_list.append(title2)

uq = list(set(list(Cell_type_list)))

Sample_Name_list = list(aux['Sample Name'])

N = len(Sample_Name_list)


#construct the gene list
gene_list = []
sample = Sample_Name_list[0]
found = False
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
    if idx_raw > 0: #first line is 
        raw_data = raw_data.split()[0]  #first column is the gene name
        gene_list.append(raw_data.decode('utf-8'))

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
            if gene == gene_list[idx_gene]:
                MATRIX[idx_gene, idx_cell] = float(raw_data[1].decode('utf-8'))
            else:
                print('Gene Does not match')



#################
#################
#Write data for full data
#cell_type_unique = list(set(list(cell_type_full_list)))  #used if the cell_type makes sense
print('Writing the Full Data>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
cell_type_unique = uq

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
        
        
        

  
