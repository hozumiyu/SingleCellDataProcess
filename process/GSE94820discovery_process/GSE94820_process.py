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


data = 'GSE94820'

inpath = "./RAW/"
outpath = "../../%s/"%data
try:
    os.makedirs(outpath)
except:
    print()

#list of all the files

#meta data
aux = pd.read_csv('GSE94820_sample.csv')
raw_data_file = open(inpath + 'GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt')
raw_data_lines = raw_data_file.readlines()
raw_data_file.close()


#ge the exp name from the raw data, and the cell type and sample name from aux
Exp_Name_list = raw_data_lines[0].split()
Cell_type_list = []
Sample_Name_list = []
#map the exp name to class names
for idx, exp in enumerate(Exp_Name_list):
    for idx_title, title in enumerate(aux['Title']):
        if exp == title:
            found = True
            break
    cell_type =  title.split('_')
    cell_type = cell_type[0]
    
    Cell_type_list.append(cell_type)
    Sample_Name_list.append(aux['Accession'][idx_title])
    if found == False:
        print('WARNING', exp, 'not found in aux!')
        break

if len(Sample_Name_list) != len(list(set(Sample_Name_list))):
    print('WARNING: There is dupolicates in the sample name. consider reprocessing')
    
#construct the gene list
gene_list = []   #the first col contains the gene info
for idx_raw, raw_data in enumerate(raw_data_lines):
    if idx_raw > 0: #first row is header
        gene = raw_data.split()[0]
        gene_list.append(gene.upper())

unique = list(set(Cell_type_list))
unique.sort()

M = len(gene_list)
N = len(Sample_Name_list)
print('Number of genes:', M)
print('Number of cells:', N)
print('Number of cell type:', len(unique))  #should be 4

if M != len(raw_data_lines) - 1:
    print('Gene length and raw data does not match')

MATRIX = np.zeros([M,N])

for idx_raw, raw_data in enumerate(raw_data_lines):
    if idx_raw > 0: #first line is the header
        idx_row = idx_raw -1 
        raw_data = raw_data.split('\t')[1:]
        
        rowdata = np.array(raw_data).astype(float)
        MATRIX[idx_row, : ] = rowdata

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
        
        

