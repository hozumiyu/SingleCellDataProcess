# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 10:50:23 2022

@author: yutah
"""

import gzip
import csv
import pandas as pd
import numpy as np

data = 'GSE59114'

label = True
gene = True
X = True

inpath = "./RAW/"
outpath = "../../%s/"%data

#PLEASE COVERT YOUR FILE TO CSV

if label:
    
    aux = pd.read_csv(data + '_sample.csv')
    Title_list = aux['Title']
    Accession_list = aux['Accession']
    
    file = open(inpath  + 'GSE59114_DBA_GEO_all.csv', 'r')
    raw_data_lines = file.readlines()
    file.close()
    
    samples = raw_data_lines[1]
    Experiment_name_list = samples.split(',')[1:]  #second line is the data, first column is just the header
        
    unique_labels = ['young_DBA_population',
             'young_DBA_LTHSC',
             'young_DBA_STHSC',
             'young_DBA_MPP',
             'old_DBA_population',
             'old_DBA_LTHSC',
             'old_DBA_STHSC',
             'old_DBA_MPP']
    map_labels = {i: unique_labels[i] for i in range(len(unique_labels))}
    map_labels_rev = {unique_labels[i]: i for i in range(len(unique_labels))}
    file = open(outpath + data + '_labeldict.txt', 'w')
    file.writelines(str(map_labels) + '\n')
    file.writelines(str(map_labels_rev) )
    file.close()
    
    cell_count = {unique_labels[i]: 0 for i in range(len(unique_labels))}
    
    Sample_Name_list = []
    Cell_type_list = []
    Label_list = []
    for idx_e, exp in enumerate(Experiment_name_list):
        label = exp.split('_')
        label = '_'.join(label[:3])
        Cell_type_list.append(label)
        Label_list.append( map_labels_rev[label])
        cell_count[label] += 1
        found = True
        for idx_a, Accession in enumerate(Accession_list):
            if Title_list[idx_a] == exp:
                found = True
                break
        if found == False:
            print('Experiment name and accession name match not found')
        Sample_Name_list.append(Accession)
    
    print('Number of Samples:', len(Sample_Name_list), cell_count)
    
    with open(outpath + data + '_labels.csv', "w", newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['Sample Name',	'Cell type','Label', 'Experiment name'])
        for idx in range( len(Experiment_name_list)):
            writer.writerow([Sample_Name_list[idx] ,Cell_type_list[idx], Label_list[idx], Experiment_name_list[idx] ])
            
if gene:
    file = open(inpath + data + '_DBA_GEO_all.csv', 'r')
    lines = file.readlines()
    file.close()
    gene_list = [] 
    for idx, line in enumerate(lines):
        if idx > 1:  #first line is indexing
            line = line.split(',')
            gene_list.append(line[0].upper())
    with open(outpath + data + '_gene.csv', "w", newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(['Index', 'Gene'])
        for idx, gene in enumerate(gene_list):
            writer.writerow([idx, gene])
    print('Number of gene:', len(gene_list))
    

if X:
    exp_names = pd.read_csv(outpath + data + '_labels.csv')['Experiment name']
    Sample_name_list = list(pd.read_csv(outpath + data + '_labels.csv')['Sample Name'])
    genes_list = pd.read_csv(outpath + data + '_gene.csv')
    file = open(inpath + data + '_DBA_GEO_all.csv', 'r')
    raw_data_lines = file.readlines()
    file.close()
    sample_names = lines[1].split(',')
    sample_names[-1] = sample_names[-1].strip()
    
    Matrix = np.zeros([len(genes_list), len(exp_names)])
    for idx_col in range(len(exp_names)):
        if sample_names[idx_col+1] != exp_names[idx_col]:
            print('Names do not match', sample_names[idx_col], exp_names[idx_col])
        for idx_row in range(2,len(lines)):
            val = lines[idx_row].split(',')
            val = val[idx_col + 1].strip()
            val = float(val)
            Matrix[idx_row-2, idx_col] = val
            
            
    with open(outpath + data + '_data.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Row', 'Col', 'Val'])
        for idx_col in range(len(exp_names)):
            for idx_row in range(len(genes_list)):
                val = Matrix[idx_row, idx_col]
                if val > 0 :
                    writer.writerow([idx_row, idx_col, val])
    
    with open(outpath + data + '_X.csv', "w", newline = '')  as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow([None] + Sample_name_list)
        for idx_row in range(len(gene_list)):
            row = list(Matrix[idx_row, :])
            writer.writerow([gene_list[idx_row] ] + row)
               