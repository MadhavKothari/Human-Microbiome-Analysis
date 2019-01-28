# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 13:20:10 2018

@author: natur
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

map_v13 = 'v13_map_uniquebyPSN.txt'
map_v35 = 'v35_map_uniquebyPSN.txt'
otu_v13 = 'v13_psn_otu.genus.fixed.txt'
otu_v35 = 'v35_psn_otu.genus.fixed.txt'

def extract_location(file):
	''' Gets all the sample sites from the data set'''
	sites = []
	with open (file,'r') as f:
		for line in f:
			if line[0] != '#':
				if line.split('\t')[5] not in sites:
					sites.append(line.split('\t')[5])
	return sites

def extract_ID(file,sample):
	''' This Function is responsible for extracting the sample ID's of the user defined location'''
	id_list = []
	with open (file,'r') as f:
		for line in f:
			if line.split()[5] == sample:
				id_list.append(line.split()[0])
	return id_list

def get_all_ids(file,sites):
	'''Gets the sample IDs for all sample sites in a dictionary'''
		sample_lists = {}
		for s in sites:
			id_list = extract_ID(file,s)
			sample_lists[s] = id_list
		return sample_lists


def extract_taxa(file):
	''' Extracts the specific level of taxa of intrest and returns a list of taxa'''
	tax = []
	with open (file,'r') as f:
		for line in f:
			if line[0] != '#':
				l = line.split('\t')
				tax.append(l[-1].split(';')[-1].strip())
	return tax

def extract_colm_number(id_list,file):
	''' Function returns the coloumn numbers of our sample of intrest'''
	colm_number = []
	with open (file,'r') as f:
		line = f.readline()
		line = f.readline().split('\t')
		for i in range(len(line)):
			if line[i] in id_list:
				colm_number.append(i)
	return colm_number

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

def count_taxa_colm(matrix,tax):
    count = {}
    threshold = 0
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

def count_sample_diversity(file,sample,tax):
	'''Counts the total number of samples that were observed in each taxa that was found in that sample site'''
	count = {}
	threshold = 0
	co = matrix > threshold
	for i in range(len(co)):
	    if co[i]:
	        if tax[i] not in count:
	            count[tax[i]] = 0 
	        count[tax[i]] += matrix[i].sum()

def normalize(matrix):
	''' normalizes the given matrix along its coloums'''
	su = matrix.sum(axis = 0) > 0 
	norm_matrix = matrix[:,su]/matrix[:,su].sum(axis = 0)
	return norm_matrix

def mla(norm_matrix,tax):
	'''Calculates the mean log abundances'''
	mla = {}
	pseudo_count = 0.1/np.mean(norm_matrix.sum(axis=0))
	mean_matrix = np.mean(np.log(norm_matrix+pseudo_count), axis = 1)
	for i in range(len(tax)):
		if tax[i] not in mla:
			mla[tax[i]] = []
		mla[tax[i]].append(mean_matrix[i])
	for key in mla:
		mla[key] = np.mean(mla[key])
	dele = []
	for key in mla:
	    if key[0] != 'g':
	        dele.append(key)
	for key in dele:
	    del mla[key]
	return mla

sites = extract_location(map_v35)
sample_lists = get_all_ids(map_v35,sites)
tax_g_v35 = extract_taxa(otu_v35)

mla_35 = {}
sample_diversity = {}
for key in sites:
	print(key)
	colm_number = extract_colm_number(sample_lists[key],otu_v35)
	matrix = exteact_smpls(colm_number,otu_v35)
	norm_matrix = normalize(matrix)
	mla_35[key] = mla(norm_matrix,tax_g_v35)

def sort_combine(sorted_count,count13,count35):
	'''Combines the filtered data from two similar lists, sorted accoriding to one of the lists''' 
	final = []
	for i in reversed(sorted_count):
		if i in count13 and i in count35:
			final.append([count13[i],count35[i]])
	return final
 
heat = np.zeros((18,18)) # Creating a heatmap comparing all the different sites to every other site
for i in range(len(mla_35)):
    for j in range(len(mla_35)):
        a = sites[i]
        b = sites[j]
        if a == b:
            continue
        sorted_count = sorted(mla_35[a], key=mla_35[a].get)
        final = sort_combine(sorted_count,mla_35[a],mla_35[b])
        final_matrix = np.array(final)
        r, p_val = pearsonr(final_matrix[:,0],final_matrix[:,1]) # Pearson coefficient is used for this comparision
        heat[i,j] = r

import seaborn as sns
mask = np.zeros_like(heat)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    ax = sns.heatmap(heat, mask=mask, vmax=1, square=True,xticklabels = sites,yticklabels=sites,linewidths = 0.1)
    plt.title('Heatmap showing correlation coefficients\n between the mean log abundances at each site')
    plt.show()

'''
The next set of code calculated the mean diversity at each site, 
from which the most and least diverse sites were used to make a rank abundance plot
'''
    
combine_v35 = {}
for key in sample_lists:
	combine = {}
	print (key)
	colm_no = extract_colm_number(sample_lists[key], otu_v35)
	matrix = exteact_smpls(colm_no,otu_v35)
	norm_matrix = normalize(matrix)
	for i in range(np.shape(norm_matrix)[1]):
		count = count_taxa_colm(norm_matrix[:,i],tax_g_v35)
		combine[sample_lists[key][i]] = count
		i+=1
	combine_v35[key] = combine

diversity = {}
for key in combine_v35:
	temp = []
	for smp_id in combine_v35[key]:
		temp.append(len(combine_v35[key][smp_id]))
	diversity[key] = temp

site_mean = {} # This stored the average diversity at each site
for key in diversity:
	site_mean[key] = np.mean(diversity[key])
'''
Rank abundance plot for most and least diverse sites (diversity was calculated beased on how many different otu's were present at each site)
'''
def count_taxa(matrix,tax):
	count = {}
	threshold = -1
	co = matrix.sum(axis = 1) > threshold
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


id_list35 = extract_ID(map_v35,'Subgingival_plaque')
colm_number35 = extract_colm_number(id_list35,otu_v35)
matrix35 = exteact_smpls(colm_number35,otu_v35)
norm_matrix35 = normalize(matrix35)
count35 = count_taxa(norm_matrix35,tax_g_v35)

sorted_count35 = sorted(count35, key=count35.get)
final = sort_combine(sorted_count35,count35,count35)

final_matrix = np.array(final)

plt.loglog(final_matrix[:,1],'r.',label = 'v_35')
plt.xlabel('Genus')
plt.ylabel('Frequency')
plt.title('Rank abundance plot of v35 regions for Subgingival_plaque at Genus level')
plt.legend()

plt.savefig('RAP v_35 Subgingival_plaque genus.png')



