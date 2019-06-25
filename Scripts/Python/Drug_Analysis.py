
import csv
from pprint import pprint
import numpy as np
import pandas as pd
from pandas import DataFrame
import itertools

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

file_drugs = 'Drug Analysis.csv'
df = pd.read_csv(file_drugs, sep='\t')
# df = df[df['number'] != 0]

# set width of bar
barWidth = 0.25

#ordering the bars
df = df.sort_values(['Degree 2Neighbors Network-proximity to disease genes', 'Degree 2Neighbors Network-proximity to disease genes&without outliers'], ascending=False )
# set height of bar
bars1 = df['Degree 2Neighbors Network-proximity to disease genes'].tolist()
# print(type(bars1))
# print(bars1)
bars2 = df['Degree 2Neighbors Network-proximity to disease genes&without outliers'].tolist()
# bars3 = []

# Set position of bar on X axis
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
# r3 = [x + barWidth for x in r2]

# Make the plot
plt.bar(r1, bars1, color='#7f6d5f', width=barWidth, edgecolor='white', label='Number of drug targets (in subnetwork based on proximity to disease genes)')
plt.bar(r2, bars2, color='#557f2d', width=barWidth, edgecolor='white', label='Number of targets for every drug(in highly connected network based on ClusterONE algorithm)')
# plt.bar(r3, bars3, color='#2d7f5e', width=barWidth, edgecolor='white', label='var3')

# Add xticks on the middle of the group bars
plt.xlabel('Drugs', fontweight='bold')

def cocatenate_origin_and_name(drug_info):
    return drug_info['type'] + ' - ' + drug_info['drugs']

df['origin'] = df.apply(cocatenate_origin_and_name, axis=1)
xlabels = df['origin'].tolist()
# plt.xticks([r + barWidth for r in range(len(bars1))], xlabels, rotation= -90, fontsize=6)
plt.xticks(r1, xlabels, rotation= -90, fontsize=6)

###############
# Make the plot
bars3 = df['Nº Clusters'].tolist()
markerline, stemlines, baseline = plt.stem(r1, bars3, linefmt='grey', markerfmt='D', label='Number of clusters targeted for every drug(in subnetwork based on proximity to disease genes)')
markerline.set_markerfacecolor('#7f6d5f')
bars4 = df['Nº Clusters2'].tolist()
markerline2, stemlines2, baseline2 = plt.stem(r2, bars4, linefmt='grey', markerfmt='D', label='Number of clusters targeted for every drug(in highly connected network)')
markerline2.set_markerfacecolor('#557f2d')

# Create legend & Show graphic
plt.legend()
plt.savefig('Drug_Analysis_Barplot.svg')
#plt.show()
###############

plt.figure()
file_drugs1 = 'number_targets_of_drugs_subNetwork.csv'
df1 = pd.read_csv(file_drugs1)

bars1 = df1.groupby('Outdegree')['drugs'].count()
# bars1 = df.groupby('Degree 2Neighbors Network-proximity to disease genes')['drugs'].count()

# Set position of bar on X axis
r1 = np.arange(len(bars1))

# Make the plot
plt.bar(r1, bars1, color='#7f6d5f', width=barWidth, edgecolor='white', label='Number of drugs with N targets(in subnetwork based on proximity to disease genes)')

# Add xticks on the middle of the group bars
plt.xlabel('Number of Targets', fontweight='bold')
plt.ylabel('Number of Drugs', fontweight='bold')

xlabels = bars1.keys()
# plt.xticks([r + barWidth for r in range(len(bars1))], xlabels, rotation= -90, fontsize=6)
plt.xticks(np.arange(len(bars1)), xlabels, fontsize=10)

plt.savefig('Number of drugs with N targets_Barplot_SubNetwork.svg')
plt.figure()

file_drugs2 = 'number_targets_of_drugs_ClusteredNetwork.csv'
df2 = pd.read_csv(file_drugs2)

bars2 = df2.groupby('Outdegree')['drugs'].count()
# bars2 = df.groupby('Degree 2Neighbors Network-proximity to disease genes&without outliers')['drugs'].count()
# bars2.pop(0)

# Set position of bar on X axis
r2 = np.arange(len(bars2))

# Make the plot
plt.bar(r2, bars2, color='#557f2d', width=barWidth, edgecolor='white', label='Number of drugs with N targets(in highly connected network based on ClusterONE algorithm)')

# Add xticks on the middle of the group bars
plt.xlabel('Number of Targets', fontweight='bold')
plt.ylabel('Number of Drugs', fontweight='bold')

xlabels = bars2.keys()
# plt.xticks([r + barWidth for r in range(len(bars1))], xlabels, rotation= -90, fontsize=6)
plt.xticks(np.arange(len(bars2)), xlabels, fontsize=10)

plt.savefig('Number of drugs with N targets_Barplot_ClusteredNetwork.svg')

plt.figure()

file_drugs3 = 'number_targets_of_drugs_2Neighbor-Network.csv'
df3 = pd.read_csv(file_drugs3, sep='\t')

bars3 = df3.groupby('number')['drugs'].count()
# bars2 = df.groupby('Degree 2Neighbors Network-proximity to disease genes&without outliers')['drugs'].count()
# bars2.pop(0)

# Set position of bar on X axis
r3 = np.arange(len(bars3))

# Make the plot
plt.bar(r3, bars3, color='#2d7f5e', width=barWidth, edgecolor='white', label='Number of drugs with N targets(in network based on 2 common neighbors strategy)')

# Add xticks on the middle of the group bars
plt.xlabel('Number of Targets', fontweight='bold')
plt.ylabel('Number of Drugs', fontweight='bold')

xlabels = bars3.keys()
# plt.xticks([r + barWidth for r in range(len(bars1))], xlabels, rotation= -90, fontsize=6)
plt.xticks(np.arange(len(bars3)), xlabels, fontsize=8)

plt.savefig('Number of drugs with N targets_Barplot_2Neighbor-Network.svg')
plt.show()
