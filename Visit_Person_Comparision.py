# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 14:01:23 2018

@author: natur
"""

import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.stats.stats import pearsonr

map_v13 = 'v13_map_uniquebyPSN.txt'
map_v35 = 'v35_map_uniquebyPSN.txt'
otu_v13 = 'v13_psn_otu.genus.fixed.txt'
otu_v35 = 'v35_psn_otu.genus.fixed.txt'

'''
This code produces the graphs comparing stool samples (Can be changed) from two random people and also the same person at two different times
'''

def extract_ID(file,sample):
	''' This Function is responsible for extracting the sample ID's of the user defined location'''
	id_list = []
	with open (file,'r') as f:
		for line in f:
			if line.split()[5] == sample and line.split()[2] != 2:
				id_list.append(line.split()[0])
	return id_list

def extract_colm_number(id_list,file):
	''' Function returns the coloumn numbers of our sample of intrest'''
	colm_number = []
	id_list_final = []
	with open (file,'r') as f:
		line = f.readline()
		line = f.readline().split('\t')
		for i in range(len(line)):
			if line[i] in id_list:
				colm_number.append(i)
				id_list_final.append(line[i])
	return colm_number, id_list_final

def extract_taxa(file):
	''' Extracts the specific level of taxa of intrest and returns a list of taxa'''
	tax = []
	with open (file,'r') as f:
		for line in f:
			if line[0] != '#':
				l = line.split('\t')
				tax.append(l[-1].split(';')[-1].strip())
	return tax

def exteact_smpls(colm_number,file):
	''' This function will return a matrix of the counts of each sample of intrest'''
	matrix = []
	threshold = 10
	with open (file,'r') as f:
		for line in f:
			if line[0] != '#':
				l = line.split('\t')
				value = []
				for i in colm_number:
					if int(l[i]) < threshold:
						value.append(0)
					else:
						value.append(int(l[i]))
				matrix.append(value)
	return np.array(matrix)

def normalize(matrix):
	''' normalizes the given matrix along its coloums'''
	norm_matrix = matrix/matrix.sum(axis = 0)
	return norm_matrix

def count_taxa_colm(matrix,tax):
    count = {}
    threshold = -1
    co = matrix > threshold
    for i in range(len(co)):
        if co[i]:
            if tax[i] not in count:
                count[tax[i]] = 0 
            count[tax[i]] += matrix[i].sum()
    dele = []
    for key in count:
        if key[0] != 'g':
            dele.append(key)
    for key in dele:
        del count[key]
    return count

def sort_combine(sorted_count,count13,count35):
	'''Combines the filtered data from both the v13 and v35 regions and sorted accoriding to the v13 data''' 
	final = []
	for i in reversed(sorted_count):
		if i in count13 and i in count35:
			final.append([count13[i],count35[i]])
	return final

stool_id = extract_ID(map_v13, "Stool") # Stool can be changed to any other sample site as the user wants
stool_colm,id_list_final = extract_colm_number(stool_id, otu_v13)
stool_matrix = exteact_smpls(stool_colm,otu_v13)
norm_matrix = normalize(stool_matrix)

stool_combinations = []
for a, b in itertools.combinations(id_list_final, 2):
    stool_combinations.append([a,b])

tax_g = extract_taxa(otu_v13)

combine = {}
for i in range(len(id_list_final)):
    count = count_taxa_colm(norm_matrix[:,i],tax_g)
    combine[id_list_final[i]] = count

pseudo_count = 0.1/np.mean(stool_matrix.sum(axis=0))


li = []
for a,b in stool_combinations:
    sorted_count = sorted(combine[a], key=combine[a].get)
    final = sort_combine(sorted_count,combine[a],combine[b])
    final_matrix = np.array(np.log(final+pseudo_count))
    r, p_val = pearsonr(final_matrix[:,0],final_matrix[:,1])
    li.append(r)

def extract_rsid(file,sample):
	'''Extracts the rsid from the table to link the data to a specif person'''
	rsid = []
	with open (file,'r') as f:
		for line in f:
			if line.split()[5] == sample and line.split()[2] == '2':
				rsid.append(line.split()[1])
	return rsid

def extract_both_visit(file,sample,rsid):
	''' Using the rsid values the data of people who have been sampled twice is extracted so it can be compared'''
	visit = {} # Saves the list of pople who have had a secon visit
	allid = [] # Saves the list of people who have visited for the first time
	with open (file,'r') as f:
		for line in f:
			for r in rsid:
				if line.split()[5] == sample and line.split()[2] == '1' and line.split()[1] == r:
					if r in visit:
						visit[r].append(line.split()[0])
					else:
						visit[r] = [line.split()[0]]
					allid.append(line.split()[0])
				if line.split()[5] == sample and line.split()[2] == '2' and line.split()[1] == r:
					if r in visit:
						visit[r].append(line.split()[0])
					else:
						visit[r] = [line.split()[0]] 
					allid.append(line.split()[0])
	return visit,allid

rsid = extract_rsid(map_v13,"Stool")
visit,allid = extract_both_visit(map_v13,"Stool",rsid)
allclm,toss = extract_colm_number(allid,otu_v13)
visit_stool_matrix = exteact_smpls(allclm,otu_v13)
visit_norm_matrix = normalize(visit_stool_matrix)
v_combine = {}
for i in range(len(toss)):
    count = count_taxa_colm(visit_norm_matrix[:,i],tax_g)
    v_combine[toss[i]] = count
same_person = []
for key in visit:
	if len(visit[key]) == 2:
		same_person.append(visit[key])

vli = []
for a,b in same_person:
	try: # A try and except is used as in some cases there is data for a 2nd vist but no 1st visit data was available
	    vsorted_count = sorted(v_combine[a], key=v_combine[a].get)
	    vfinal = sort_combine(vsorted_count,v_combine[a],v_combine[b])
	    vfinal_matrix = np.array(np.log(vfinal+pseudo_count))
	    r, p_val = pearsonr(vfinal_matrix[:,0],vfinal_matrix[:,1])
	    vli.append(r)
	except:
		pass






