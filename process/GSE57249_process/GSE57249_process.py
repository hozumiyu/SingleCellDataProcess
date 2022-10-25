# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 10:50:23 2022

@author: yutah
"""

import gzip
import csv
import pandas as pd
import numpy as np
data = 'GSE57249'


label = True
gene = True
X = True
inpath = "./RAW/"
outpath = "../../%s/"%data



if label:
    meta = pd.read_csv(outpath + data + '_meta.csv')
    GEO_Accession = meta['GEO_Accession (exp)']
    
    file = open(inpath + data + '_fpkm.txt')
    lines = file.readlines()
    file.close()
    sample_name_list = lines[0].split()[1:]
    
    cell_type_list = meta['Developmental_stage']
    cell_type_unique = ['zygote', 'Two-cell Embryo Blastomere', 
                        'Four-cell Embryo Blastomere', 
                        'inner cell mass from blastocyst', 
                        'trophectoderm from blastocyst']
    cell_type_dict = {i: cell_type_unique[i] for i in range(len(cell_type_unique)) }
    cell_type_dict_rev = {cell_type_unique[i]:i for i in range(len(cell_type_unique))}
    file = open(outpath + data + '_labeldict.txt', 'w')
    file.writelines(str(cell_type_dict) + ' \n')
    file.writelines(str(cell_type_dict_rev) + ' \n')
    file.close()
    cell_count = {cell_type_unique[i]:0 for i in range(len(cell_type_unique))}
    Sample_Name_list = []
    Cell_type_list = []
    Cell_label_list = []
    for sample_idx, sample in enumerate(sample_name_list):
        found = False
        for idx, geo in enumerate(GEO_Accession):   #make sure that sample name and GEO ACC matches
            if geo == sample:
                found = True
                
        if found == False:
            print('Sample Name does not match')
        cell_type = cell_type_list[idx]
        if cell_type == '"inner cell mass\\_ from blastocyst"':  #adjust the name of the cells
            cell_type = 'inner cell mass from blastocyst'
        elif cell_type == '"trophectoderm\\_ from blastocyst"':
            cell_type = 'trophectoderm from blastocyst'
        cell_count[cell_type] += 1
        Cell_type_list.append(cell_type); Cell_label_list.append(cell_type_dict_rev[cell_type]); Sample_Name_list.append(GEO_Accession[idx])
    with open(outpath + data + '_labels.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sample Name', 'Cell type', 'Label'])
        for idx in range(len(Sample_Name_list)):
            writer.writerow([Sample_Name_list[idx], Cell_type_list[idx], Cell_label_list[idx]])
    print('Number of Samples:', len(Sample_Name_list))
    print(cell_count)

if gene:
    file = open(inpath + data + '_fpkm.txt')
    lines = file.readlines()
    file.close()
    gene_list = [] 
    for idx, line in enumerate(lines):
        if idx > 0:  #first line is indexing
            line = line.split()
            gene_list.append(line[0])
    with open(outpath + data + '_gene.csv', "w", newline = '') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['Index', 'Gene'])
        for idx, gene in enumerate(gene_list):
            writer.writerow([idx, gene])
    print('Number of gene:', len(gene_list))
    
            
if X:
    Sample_Name_list = list(   pd.read_csv(outpath + data + '_labels.csv')['Sample Name']   )
    gene_list = list(  pd.read_csv(outpath + data + '_gene.csv')['Gene']  )
    
    raw_data_file = open(inpath + data + '_fpkm.txt', 'r')
    raw_data_lines = raw_data_file.readlines()
    raw_data_file.close()
    
    Matrix = np.zeros([len(gene_list), len(Sample_Name_list)])
    for idx_raw, raw_data in enumerate(raw_data_lines):
        idx_row = idx_raw - 1
        if idx_raw > 0:
            raw_data = raw_data.split()[1:]
            for idx_col, val in enumerate(raw_data):
                val = float(val)
                if val > 0:
                    Matrix[idx_row, idx_col] = val
        
    with open(outpath + data + '_data.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Row', 'Col', 'Val'])
        for idx_col in range(len(Sample_Name_list)):
            for idx_row in range(len(gene_list)):
                val = Matrix[idx_row, idx_col]
                if val > 0 :
                    writer.writerow([idx_row, idx_col, val])
    with open(outpath + data + '_X.csv', "w", newline = '') as csv_file: #output file for data
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow( [None] + Sample_Name_list)
        for idx_row in range(len(gene_list)):
            row = list(Matrix[idx_row, :])
            writer.writerow([gene_list[idx_row] ] + row)
               
