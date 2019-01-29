# Human-Microbiome-Analysis
Codes used to extract data and perfom further analysis from the human microbiome project

The codes all share some common functions, the second half of the code is used to call these functions in the required order to extract the data and plot it as required. The input files are also inculded here. They were obtained from the human microbiome project site (https://hmpdacc.org/hmp/HMQCP/). The files numbered 9 and 10 were taken for this analysis. 

Each code focuses on a specifc aspect.

MLA and RAP:
-This code compares the mean log abundances and the rank abundance plot between the data obtasined from the v13 and v35 regions. This analysis was run on different OTU levels (family and genus, varaible can be changed for other OTU levels). 

Visit Person Comparision:
-This runs a comparision between two random people on the same visit, and also a compariosion between a random person between his/her two visits. The code also produces an histogram of the pearson coefficient between all the samples of the same visit, and a histogram of the perason coefficient of the two visits of the same person.  

Site Diversity and Site Comparision:
-This compares all the sites to one another to produce a heatmap of how similar each site is to every other site. This similarity is based on the perarson coeffcient of the rank abundance of the site. The comparision is done by sorting both the sites according to only one of the sites. This also produces a table that has the calculated mean diversity of each site. The most and least diverse of all the sites were then used to produce a rank abundance plot.  

The img folder contains some of the types of graphs that these codes will generate. 
