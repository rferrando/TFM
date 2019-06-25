
import csv
from pprint import pprint
import numpy as np
import pandas as pd
from pandas import DataFrame
import itertools

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# file_drugs1 = 'number_targets_of_drugs_subNetwork.csv'
file_drugs1 = 'number_targets_of_drugs_2Neighbor-Network.csv'
df = pd.read_csv(file_drugs1, sep='\t')
# print(df.dtypes)
##################################

print("Drug Boxplot")
# print(df.dtypes)
data_plot = df['Outdegree']
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
ax.set_title('Distribución de drogas por nº de interacciones')
ax.set_xticklabels(['Drogas cuyos target son genes en la subred generada en base a la proximidad a genes de enfermedad'])
ax.set_ylabel('Nº de targets')

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
plt.savefig('Drug_Boxplot_2Neighbor-Network.svg')
# plt.figure()

print("Drug Boxplot by Type")

data_plot = [df[df['type'] == 'E']['Outdegree'], df[df['type'] == 'I']['Outdegree'], df[df['type'] == 'S']['Outdegree']]
fig, ax = plt.subplots(figsize=(20, 10))
bp = ax.boxplot(data_plot)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')
# colors = ['darkkhaki', 'royalblue', 'salmon']
colors = ['lightgreen', 'tan', 'pink',]

ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

# Hide these grid behind plot objects
ax.set_axisbelow(True)
ax.set_title('Comparación del nº de targets para cada tipo de droga')
ax.set_xlabel('Tipo de droga')
ax.set_ylabel('Nº de targets')
num_boxes = len(data_plot)
for i in range(num_boxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    box_coords = np.column_stack([boxX, boxY])
    # Alternate between Dark Khaki and Royal Blue
    ax.add_patch(Polygon(box_coords, facecolor=colors[i % 3]))
# Multiple box plots on one Axes
ax.set_xticklabels(['Drogas cuyos targets son genes causantes de endometriosis','Drogas que interactuan con drogas de endometriosis', 'Drogas relacionadas con endometriosis'],
                    fontsize=9)#rotation=45, fontsize=8)
# Calculate number of obs per group & median to position labels
medians = df.groupby(['type'])['Outdegree'].median().sort_index(axis=0).values
# print('medians: ', df.groupby(['type'])['Outdegree'].median().sort_index(axis=0).keys())
# medians:  Index(['E', 'I', 'S'], dtype='object', name='type')
print(medians)
nobs = df['type'].value_counts().sort_index(axis=0).values
# print('nobs: ', df['type'].value_counts().sort_index(axis=0).keys())
# nobs:  Index(['E', 'I', 'S'], dtype='object')


nobs = [str(x) for x in nobs.tolist()]
nobs = ["n: " + i for i in nobs]
# nobs.reverse()
# Add it to the plot
pos = range(len(nobs))
for tick,label in zip(pos,ax.get_xticklabels()):
    # print(tick, label)
    # 0 Text(0, 0, 'Drogas cuyos targets son genes causantes de endometriosis') -> 'E'
    # 1 Text(0, 0, 'Drogas que interactuan con drogas de endometriosis') -> 'I'
    # 2 Text(0, 0, 'Drogas relacionadas con endometriosis') -> 'S'

    ax.text(pos[tick] + 1, medians[tick] - 0.25, nobs[tick],
    horizontalalignment='center', size='x-small', color='black', weight='semibold')

# plt.savefig('Drug_Boxplot_subnetwork_by_type.svg')
plt.savefig('Drug_Boxplot_2Neighbor-Network_by_type.svg')

print("Drug Density Plot")
import seaborn as sns
# plot of 2 variables
# plt.figure()
p1=sns.kdeplot(df[df['type'] == 'S']['Outdegree'], shade=True, color=colors[2], label ='Drogas relacionadas con endometriosis')
p2=sns.kdeplot(df[df['type'] == 'I']['Outdegree'], shade=True, color=colors[1],label ='Drogas que interactuan con drogas de endometriosis')
p3=sns.distplot(df[df['type'] == 'E']['Outdegree'], color=colors[0], label ='Drogas cuyos targets son genes causantes de endometriosis', hist=False)
# plt.savefig('Drug_Density_plot_Subnetwork.svg')
plt.savefig('Drug_Density_plot_2Neighbor-Network.svg')
plt.show()
