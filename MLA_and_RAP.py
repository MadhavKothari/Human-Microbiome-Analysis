# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 15:38:14 2018

@author: Naturedrag
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

'''
This code is used to make a mean log abundance plot and a rank abundance plot 
from one sample site in both the v13 and v35 regions to compare them. The plots 
are made at the genus level, but can be changed.
'''

map_v13 = 'v13_map_uniquebyPSN.txt'
map_v35 = 'v35_map_uniquebyPSN.txt'
otu_v13 = 'v13_psn_otu.genus.fixed.txt'
otu_v35 = 'v35_psn_otu.genus.fixed.txt'

def extract_ID(file,sample):
	''' This Function is responsible for extracting the sample ID's of the user defined location'''
	id_list = []
	with open (file,'r') as f:
		for line in f:
			if line.split()[5] == sample and line.split()[2] != '2':
				id_list.append(line.split()[0])
	return id_list

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

#def extract_taxa(file): # This function is used to calculte for any taxanomic level other than genus 
#	''' Extracts the specific level of taxa of intrest and returns a list of taxa '''
#	tax = []
#	with open (file,'r') as f:
#		for line in f:
#			if line[0] != '#':
#				l = line.split('\t')
#				if len(l[-1].split(';')) > 4:
#					tax.append(l[-1].split(';')[4].strip())
#				else:
#					tax.append('N/A')
#	return tax

def extract_taxa(file):
	''' Extracts the specific level of taxa of intrest and returns a list of taxa'''
	tax = []
	with open (file,'r') as f:
		for line in f:
			if line[0] != '#':
				l = line.split('\t')
				tax.append(l[-1].split(';')[-1])
	return tax

def count_taxa(matrix,tax):
	'''Counts the number of times each taxa has appeard in the table'''
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
		if key[0] != 'g': # This has to be changed when the taxanomic level is changed
			dele.append(key)
	for key in dele:
		del count[key]
	return count

def normalize(matrix):
	''' normalizes the given matrix along its coloums'''
	norm_matrix = matrix/matrix.sum(axis = 0)
	return norm_matrix

def sort_combine(sorted_count,count13,count35):
	'''Combines the filtered data from both the v13 and v35 regions and sorted accoriding to the v13 data''' 
	final = []
	for i in reversed(sorted_count):
		if i in count13 and i in count35:
			final.append([count13[i],count35[i]])
	return final

def mla(norm_matrix,tax):
	'''Calculates the mean log abundances of the different otu's (at a user defined taxanomic level)'''
	mla = {}
	pseudo_count = 0.1/np.mean(norm_matrix.sum(axis=0))
	mean_matrix = np.mean(np.log(norm_matrix+pseudo_count), axis = 1)
	for i in range(len(tax)):
		if tax[i] not in mla:
			mla[tax[i]] = []
		mla[tax[i]].append(mean_matrix[i])
	for key in mla:
		mla[key] = np.mean(mla[key])
	return mla

# The following code segment performs the steps requred to extract the data and filter it so as to plot the required dataset				
id_list13 = extract_ID(map_v13,'Tongue_dorsum') # The Tounge_dorsum can be changed to obtain the data from the different sample sites
colm_number13 = extract_colm_number(id_list13,otu_v13)
matrix13 = exteact_smpls(colm_number13,otu_v13)
tax13 = extract_taxa(otu_v13)
norm_matrix13 = normalize(matrix13)
count13 = count_taxa(norm_matrix13,tax13)
mla13 = mla(norm_matrix13,tax13)

#This segment does the same as the above code, except for the v35 data

id_list35 = extract_ID(map_v35,'Tongue_dorsum')
colm_number35 = extract_colm_number(id_list35,otu_v35)
matrix35 = exteact_smpls(colm_number35,otu_v35)
tax35 = extract_taxa(otu_v35)
norm_matrix35 = normalize(matrix35)
count35 = count_taxa(norm_matrix35,tax35)
mla35 = mla(norm_matrix35,tax35)

#this combines both the data sets to make a final matrix
sorted_count13 = sorted(count13, key=count13.get)
final = sort_combine(sorted_count13,count13,count35)

final_matrix = np.array(final)
print (pearsonr(final_matrix[:,0],final_matrix[:,1]))

# The remainder of the code is used to plot the required graphs
plt.loglog(final_matrix[:,0],'b.',label = 'v_13')
plt.loglog(final_matrix[:,1],'r.',label = 'v_35')
plt.xlabel('Genus')
plt.ylabel('Frequency')
plt.title('Rank abundance plot of v13 vs v35 regions\n for Tounge samples at Genus level')
plt.legend()

plt.savefig('RAP v_13 vs v_35 tounge genus.png')

mlasorted_count13 = sorted(mla13, key=mla13.get)
mlafinal = np.matrix(sort_combine(mlasorted_count13,mla13,mla35))
print (pearsonr(mlafinal[:,0],mlafinal[:,1]))

plt.semilogx(mlafinal[:,0],'b.', label = 'v_13')
plt.semilogx(mlafinal[:,1],'r.', label = 'v_35')
plt.xlabel('Family')
plt.ylabel('Mean Log Abundance')
plt.title('Mean Log Abundance of v_13 vs v_35 regions\n for Toungle Samples at Family level')
plt.savefig('MLA v_13 vs v_35 tounge family.png')





