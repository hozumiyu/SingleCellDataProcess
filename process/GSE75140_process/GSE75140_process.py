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


data = 'GSE75140'

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
raw_data_file = open(inpath + 'GSE75140_hOrg.fetal.master.data.frame.txt')
raw_data_lines = raw_data_file.readlines()
raw_data_file.close()


#construct the experiment list from the raw data
Exp_Name_list = []
for idx_raw, raw_data in enumerate(raw_data_lines):
    if idx_raw > 0: #first line is the header
        raw_data = raw_data.split()
        Exp_Name_list.append(eval(raw_data[0]))

#Construct the sample list from aux file
Sample_Name_list = []
for exp in Exp_Name_list:
    found = False
    for idx_title, title in enumerate(aux['Title']):
        title = title.split('_')
        title = '_'.join(title[1:])
        if exp == title:
            found = True
            break
    if found == False and exp == 'F5_fetal_12wpc_c1':
        Sample_Name_list.append('GSM1957793')
    elif found == False:
        print('Cannot find experiment name in aux file', exp)
        break
    else:
        Sample_Name_list.append(aux['Accession'][idx_title])

for sample in meta['GEO_Accession (exp)']:
    if sample not in Sample_Name_list:
        print(sample)

Cell_type_list = []  #construct the cell list from STAGE in meta
for sample in Sample_Name_list:
    found = False
    for idx_geo, geo in enumerate(meta['GEO_Accession (exp)']):
        if geo == sample:
            found = True
            break
    if found == False:
        print('Cannot find geo in meta data')
    if sample == 'GSM1957673' and meta['STAGE'][idx_geo] == '-':
        print('Correcting cell type for', sample)
        Cell_type_list.append('12 weeks post-conception')
    else:
        Cell_type_list.append(meta['STAGE'][idx_geo])


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



unique = ['12 weeks post-conception',
            '13 weeks post-conception',
            '33 days',
            '35 days',
            '37 days',
            '41 days',
            '53 days',
            '58 days',
            '65 days']

N = len(Sample_Name_list)
print('Number of cells:', N)
print('Cell types:', len(unique))


#construct the gene list
gene_list = raw_data_lines[0]   #the first row contains the gene 
gene_list = gene_list.split()[1:-1]
for idx_gene, gene in enumerate(gene_list):
    gene_list[idx_gene] = eval(gene).upper()

print('Number of Genes:', len(gene_list))

MATRIX = np.zeros([len(gene_list), N])


for idx_raw, raw_data in enumerate(raw_data_lines):
    idx_col = idx_raw-1
    if idx_raw > 0:
        raw_data = raw_data.split()[1:-1]
        data_col = np.zeros([len(raw_data)])
        if len(raw_data) != len(gene_list):
            print('Number of entries does not match')
            break
        for idx_gene in range(len(raw_data)):
            data_col[idx_gene] = float(raw_data[idx_gene])
        MATRIX[:, idx_col] = data_col

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
        
        

  
