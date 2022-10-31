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


data = 'GSE84133'


inpath = "./RAW/"


#list of all the files
exp_list = ['human1', 'human2', 'human3', 'human4', 'mouse1', 'mouse2']
onlyfiles = [f for f in listdir(inpath) if isfile(join(inpath, f))]
#meta data
for exp in exp_list:
    for file in onlyfiles:
        if file[-2:] == 'gz':
            if exp in file:
                raw_data = pd.read_csv(inpath + file)
                outpath = '../../%s%s/'%(data, exp)
                try:
                    os.makedirs(outpath)
                except:
                    print()
                #headers are the genes
                gene_list = list(raw_data.keys())[3:]
                for idx_gene, gene in enumerate(gene_list):
                    gene_list[idx_gene] = gene.upper()
                if len(gene_list) != len(list(raw_data.keys())) - 3:
                    print('Gene list size does not match raw data')
                Sample_Name_list = list(raw_data['Unnamed: 0'])
                Cell_type_list = list(raw_data['assigned_cluster'])
                unique = list(set(Cell_type_list))
                unique.sort()
                

                N = len(Sample_Name_list)
                M = len(gene_list)
                
                if len(raw_data) != N:
                    print('Number of cells do not match')
                
                print('Current exp:', exp)
                print('Number of cells:', N)
                print('Number of Genes:', M)
                print('Number of cell type:', len(unique))
                
                MATRIX = raw_data.values[:, 3:].astype(float)
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
                file = open(outpath + data +exp+ '_full_labeldict.txt', 'w')   #write a dictionary
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
                
                with open(outpath + data + exp+ '_full_labels.csv', 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(['Sample Name', 'Cell type', 'Label'])
                    for idx in range(len(Cell_type_list)):
                        writer.writerow([Sample_Name_list[idx], Cell_type_list[idx] , Cell_label_list[idx]])
                        
                                
                with open(outpath + data + exp+'_full_gene.csv', "w", newline = '') as csv_file:
                    writer = csv.writer(csv_file, delimiter=',')
                    writer.writerow(['Index', 'Gene'])
                    for idx, gene in enumerate(gene_list):
                        writer.writerow([idx, gene])
                                
                
                
                with open(outpath + data + exp+ '_full_data.csv', "w", newline = '') as csv_file: #output file for data
                    writer = csv.writer(csv_file, delimiter=',')
                    writer.writerow(['Row', 'Col', 'Val'])
                    for idx_col in range(len(Sample_Name_list)):
                        for idx_row in range(len(gene_list)):
                            if MATRIX[idx_row, idx_col] > 0:
                                writer.writerow([idx_row, idx_col, MATRIX[idx_row, idx_col]])
                                    
                with open(outpath + data +exp+ '_full_X.csv', "w", newline = '') as csv_file: #output file for data
                    writer = csv.writer(csv_file, delimiter=',')
                    writer.writerow( [None] + Sample_Name_list)
                    for idx_row in range(len(gene_list)):
                        row = list(MATRIX[idx_row, :])
                        writer.writerow([gene_list[idx_row] ] + row)
                        
                print('>>>>>>>>>>>>>>>>>>>>>>>')
                
                  
