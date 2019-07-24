
import sys
import csv
from pprint import pprint
import numpy as np
import pandas as pd
from pandas import DataFrame
import itertools

print("Reading Annotation Variants File")

file_variants = "/home/rosa/Master_Bioinformatica/TFM/HGMD/hgmd/hgmd_pro_2019.1_hg38_annot_filter.vcf"
# file_variants = "/home/rosa/Master_Bioinformatica/TFM/HGMD/hgmd/hgmd_pro_2019.1_hg38_annot_filter_aux.vcf"

variants = []


def AddNewPair(dict, key, valor):
     dict.update({key: valor});
     return dict

def parse_info(field, dict):
    pair = field.split('=')
    dict = AddNewPair(dict, pair[0], pair[1])
    return dict

maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

def get_impact(transcript):
    return transcript.split('|')[2]


with open(file_variants, "r") as a:
    reader = csv.reader(a, delimiter='\t',doublequote='False', quotechar='', quoting=csv.QUOTE_NONE)
    test = True
    for row in reader:
        if str(row[0])[0] == '#':
            # print(row)
            continue

        info = row[7].split(';')
        attributes = {}

        for field in info:
            attributes = parse_info(field, attributes)
#Obtener el impacto de todos los transcritos
        attributes['IMPACT'] = attributes['ANN'].split(',')
        map_iterator = map(get_impact, attributes['IMPACT'])
        attributes['IMPACT'] = list(map_iterator)
        # attributes['IMPACT'] = ['MODIFIER', 'MODIFIER', 'MODIFIER', 'MODIFIER', 'MODIFIER', 'MODIFIER', 'MODIFIER']
        attributes['IMPACT'] = dict(zip(attributes['IMPACT'],map(attributes['IMPACT'].count,attributes['IMPACT'])))
# Header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        attributes = AddNewPair(attributes, 'HGMD Variant', row[2])
        attributes = AddNewPair(attributes, 'CHROM', row[0])
        attributes = AddNewPair(attributes, 'POS', row[1])
        attributes = AddNewPair(attributes, 'REF', row[3])
        attributes = AddNewPair(attributes, 'ALT', row[4])
        # if test == True:
        #     test = False
        #     print(row)

        variants.append(attributes)

# print(variants[5])
# print(variants[6])

rankscore = []
moderate = []
high = []
low = []
modifier = []
rankscore_high = []
rankscore_moderate = []
rankscore_low = []
rankscore_modifier = []

rankscore_h = []
rankscore_m = []
rankscore_l = []

def get_attribute_by_variant(attribute, value):
    if any(value in d for d in attribute):
        return float(attribute[value])
    else:
        if value == 'RANKSCORE':
            return 2
        else:
            return 0

def get_rankscore_by_aggregate_impact(impact, values):
    impact_values_set = {'HIGH', 'MODERATE', 'LOW'}
    values_set = set()
    values_set.update(values)
    other_set_impact = list(impact_values_set.difference(values_set))[0]
    # print(other_set_impact)
    # print(any(values[0] in d for d in impact))
    # print(not any(other_set_impact in d for d in impact))
    # print(any(values[1] in d for d in impact))
    if any(values[0] in d for d in impact) and not any(other_set_impact in d for d in impact):
        return True
    else:
        if len(values) > 1 and any(values[1] in d for d in impact) and not any(other_set_impact in d for d in impact):
            return True
        else:
            return False

def get_rankscore_by_impact(impact, values):
    impact_values_set = {'HIGH', 'MODERATE', 'LOW'}
    values_set = set()
    values_set.update(values)
    other_set_impact = list(impact_values_set.difference(values_set))
    # print(other_set_impact)
    # print(any(values[0] in d for d in impact))
    # print(not any(other_set_impact in d for d in impact))
    # print(any(values[1] in d for d in impact))
    if any(values[0] in d for d in impact) and not any(other_set_impact[0] in d for d in impact) and not any(other_set_impact[1] in d for d in impact):
        return True
    else:
        return False

p = [0,0,0]
number_of_variants_with_more_than_one_impact = 0
df_high_variants = pd.DataFrame()

variants_header = ['ALT', 'ANN', 'CHROM', 'CLASS', 'DNA', 'GENE', 'HGMD Variant', 'IMPACT', 'LOF', 'MUT', 'NMD', 'PHEN', 'POS', 'PROT', 'RANKSCORE', 'REF', 'STRAND', 'DB']
data_for_variants = {}
n = 0
for variant in variants:
    # print(variant)
    n = n + 1
    key = 'row' + str(n)
    data_for_variants.update({key: variant})
#### Rankscore Boxplot for every impact
    # if get_rankscore_by_aggregate_impact(variant['IMPACT'], ['HIGH', 'MODERATE']):
    #     rankscore_h.append(get_attribute_by_variant(variant, 'RANKSCORE'))
    # if get_rankscore_by_aggregate_impact(variant['IMPACT'], ['LOW','MODERATE']):
    #     rankscore_l.append(get_attribute_by_variant(variant, 'RANKSCORE'))
    with_more_than_one_impact = False
    if get_rankscore_by_impact(variant['IMPACT'], ['HIGH']):
        rk = get_attribute_by_variant(variant, 'RANKSCORE')

        if rk != 2:
            rankscore_h.append(rk)
            df_high_variants = df_high_variants.append(variant, ignore_index=True)
        else:
            p[0] = p[0] + 1
            # rankscore_h.append(rk)
            df_high_variants = df_high_variants.append(variant, ignore_index=True)
    else:
        # print(variant)

        if get_rankscore_by_impact(variant['IMPACT'], ['MODERATE']):
            rk = get_attribute_by_variant(variant, 'RANKSCORE')

            if rk != 2:
                rankscore_m.append(rk)
            else:
                p[1] = p[1] + 1
                # rankscore_m.append(rk)
        else:
            # print(variant)

            if get_rankscore_by_impact(variant['IMPACT'], ['LOW']):
                rk = get_attribute_by_variant(variant, 'RANKSCORE')

                if rk != 2:
                    rankscore_l.append(rk)
                else:
                    p[2] = p[2] + 1
                    # rankscore_l.append(rk)
            else:
                # print(variant)
                with_more_than_one_impact = True

    if (with_more_than_one_impact == True):
        number_of_variants_with_more_than_one_impact = number_of_variants_with_more_than_one_impact + 1


####    Correlation Rankscore-Impact
    #rankscore.append(get_attribute_by_variant(variant, 'RANKSCORE'))
    # moderate.append(get_attribute_by_variant(variant['IMPACT'], 'MODERATE'))
    # high.append(get_attribute_by_variant(variant['IMPACT'], 'HIGH'))
    # low.append(get_attribute_by_variant(variant['IMPACT'], 'LOW'))
    # modifier.append(get_attribute_by_variant(variant['IMPACT'], 'MODIFIER'))

#### Rankscore Probability Density

    if get_attribute_by_variant(variant['IMPACT'], 'HIGH') > 0:
        rankscore_high.append(get_attribute_by_variant(variant, 'RANKSCORE'))
    if get_attribute_by_variant(variant['IMPACT'], 'MODERATE') > 0:
        rankscore_moderate.append(get_attribute_by_variant(variant, 'RANKSCORE'))
    if get_attribute_by_variant(variant['IMPACT'], 'LOW') > 0:
        rankscore_low.append(get_attribute_by_variant(variant, 'RANKSCORE'))
    if get_attribute_by_variant(variant['IMPACT'], 'MODIFIER') > 0:
        rankscore_modifier.append(get_attribute_by_variant(variant, 'RANKSCORE'))

# print(df_high_variants.iloc[0])
file_high_variants = 'variants_with_high_impact_regardless_rankscore.csv'
df_high_variants.to_csv(file_high_variants, index=False, sep='\t')

# >>> data = {'row_1': [3, 2, 1, 0], 'row_2': ['a', 'b', 'c', 'd']}
df_variants =  pd.DataFrame.from_dict(data_for_variants, orient='index', columns=variants_header)
# print(df_variants.iloc[0])
df_variants = df_variants[df_variants['PHEN'].str.contains('endometriosis')]
# print(df_variants.shape)
file_variants_endometriosis = 'variants_endometriosis.csv'
df_variants.to_csv(file_variants_endometriosis, index=False, sep='\t')

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import seaborn as sns
import numpy as np


colors = ['darkkhaki', 'royalblue', 'salmon', 'brown']
# print("Impact Correlation Plot")
# plt.figure()
# print('Correlation Rankscore-HIGH Impact')
# print(np.corrcoef(high, rankscore))
# matplotlib.style.use('ggplot')
#
# plt.scatter (high, rankscore)
# plt.savefig('Rankscore-HIGH_Impact_Correlation_plot.svg')
#
# plt.figure()
# print('Correlation Rankscore-MODERATE Impact')
# print(np.corrcoef(moderate, rankscore))
# matplotlib.style.use('ggplot')
#
# plt.scatter (moderate, rankscore)
# plt.savefig('Rankscore-MODERATE_Impact_Correlation_plot.svg')
#
# plt.figure()
# print('Correlation Rankscore-LOW Impact')
# print(np.corrcoef(low, rankscore))
# matplotlib.style.use('ggplot')
#
# plt.scatter (low, rankscore)
# plt.savefig('Rankscore-LOW_Impact_Correlation_plot.svg')
#
# plt.figure()
# print('Correlation Rankscore-MODIFIER Impact')
# print(np.corrcoef(modifier, rankscore))
# matplotlib.style.use('ggplot')
#
# plt.scatter (modifier, rankscore)
# plt.savefig('Rankscore-MODIFIER_Impact_Correlation_plot.svg')
# plt.figure()
print("Rankscore Density Plot")

# plot of 2 variables

p1=sns.kdeplot(rankscore_high, shade=True, color=colors[0], label ='Rankscore para variantes con algún transcrito HIGH')
p1=sns.kdeplot(rankscore_moderate, shade=True, color=colors[1], label ='Rankscore para variantes con algún transcrito MODERATE')
p1=sns.kdeplot(rankscore_low, shade=True, color=colors[2], label ='Rankscore para variantes con algún transcrito LOW')
p1=sns.kdeplot(rankscore_modifier, shade=True, color=colors[3], label ='Rankscore para variantes con algún transcrito MODIFIER')

plt.savefig('Rankscore_Density_plot.svg')
#plt.figure()

print("Impact Boxplot")


data_plot = [rankscore_h, rankscore_m, rankscore_l]
# data_plot = [rankscore_h, rankscore_l]
fig, ax = plt.subplots(figsize=(20, 10))
bp = ax.boxplot(data_plot)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')

ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

# Hide these grid behind plot objects
ax.set_axisbelow(True)
ax.set_title('Comparación del rankscore por nivel de impacto')
ax.set_xlabel('Nivel de impacto')
ax.set_ylabel('Rankscore')
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
ax.set_xticklabels(['HIGH','MODERATE','LOW'],
# ax.set_xticklabels(['HIGH-MODERATE','MODERATE-LOW'],
                    fontsize=9)#rotation=45, fontsize=8)
# Calculate number of obs per group & median to position labels

import statistics
medians = [statistics.median(rankscore_h),statistics.median(rankscore_m),statistics.median(rankscore_l)]
nobs = [len(rankscore_h), len(rankscore_m), len(rankscore_l)]
# medians = [statistics.median(rankscore_h),statistics.median(rankscore_l)]
# nobs = [len(rankscore_h), len(rankscore_l)]
nobs = [str(x) for x in nobs]
nobs = ["n: " + i for i in nobs]
# nobs.reverse()
# Add it to the plot
pos = range(len(nobs))

only_one_impact_number_of_variants = 0
for tick,label in zip(pos,ax.get_xticklabels()):
    print(label, pos[tick] + 1, medians[tick], nobs[tick])

    only_one_impact_number_of_variants = only_one_impact_number_of_variants + int(nobs[tick].split('n: ')[1])
    ax.text(pos[tick] + 1, medians[tick] - 0.03, nobs[tick],
    horizontalalignment='center', size='x-small', color='w', weight='semibold')

print('Total de variantes', len(variants))
pt = p[0] + p[1] + p[2]
print('Total tratadas', pt + number_of_variants_with_more_than_one_impact + only_one_impact_number_of_variants)
print('Sin rankscore', p, ' = ', pt)
print('Con rankscore y un solo impacto por variante', only_one_impact_number_of_variants)
print('Con rankscore pero con más de un impacto por variante o MODIFIER', number_of_variants_with_more_than_one_impact)
plt.savefig('3Impact_Boxplot_without_outliers.svg')


plt.show()
