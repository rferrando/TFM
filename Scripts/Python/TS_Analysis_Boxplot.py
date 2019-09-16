
import csv
from pprint import pprint
import numpy as np
import pandas as pd
from pandas import DataFrame
import itertools

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# file_drugs1 = 'SIF_drogas.csv'
# file_drugs1 = 'SIF_final.csv'
file_drugs1 = 'SIF_TPConnections.csv'
df = pd.read_csv(file_drugs1)# sep='\t')
# print(df.dtypes)
##################################
# df = df[df['DiseaseConnections'] != 0]
# df = df[df['TPConnections'] != 0]

print("Drug Boxplot")
# print(df.dtypes)
# data_plot = df[(df['primary'] == 'None') & (df['target'] == 'True')]['DiseaseConnections']
# data_plot = df[(df['primary'] == 'None') & (df['target'] == 'True')]['TPConnections']
data_plot = df[df['disease'] == 'True']['TPConnections']
fig, ax = plt.subplots(figsize=(20, 10))


bp = ax.boxplot(data_plot, notch=True, patch_artist=True)

# colors = ['lightblue', 'lightgreen', 'tan', 'pink', 'cyan']
# for patch, color in zip(bp['boxes'], colors):
#     patch.set_facecolor(color)

plt.setp(bp['boxes'], color='lightblue')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')

ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

# Hide these grid behind plot objects
ax.set_axisbelow(True)
# ax.set_title('Distribución de Secondary Targets por nº de interacciones con Disease genes')
# ax.set_title('Distribución de Secondary Targets por nº de interacciones con Primary Targets genes')
ax.set_title('Distribución de Disease por nº de interacciones con Primary Targets genes')
# ax.set_xticklabels(['Secondary Targets en la red final de 2 vecinos'])
ax.set_xticklabels(['Disease en la red final de 2 vecinos'])
# ax.set_ylabel('Nº de Disease genes')
ax.set_ylabel('Nº de Primary Targets')

# Calculate number of obs per group & median to position labels
import statistics
medians = [statistics.median(data_plot)]
# print(medians)
nobs = [len(data_plot)]
# print(nobs)

nobs = [str(x) for x in nobs]
nobs = ["n: " + i for i in nobs]
nobs.reverse()
# Add it to the plot
pos = range(len(nobs))

for tick,label in zip(pos,ax.get_xticklabels()):
    # print(label, pos[tick] + 1, medians[tick], nobs[tick])
    ax.text(pos[tick] + 1, medians[tick] + 0.06, nobs[tick],
    horizontalalignment='center', size='x-small', color='black', weight='semibold')

# plt.savefig('Drug_Boxplot_SubNetwork.svg')
yval = np.concatenate([line.get_ydata() for line in bp['whiskers']])
print(yval.max())
# eps = 1.0
# ymin, ymax = yval.min()-eps, yval.max()+eps
# ax.set_ylim([ymin,ymax])
# plt.savefig('Boxplot_TS_connected_to_TP.svg')
plt.savefig('Boxplot_D_connected_to_TP_con0.png')
# plt.figure()
plt.show()
