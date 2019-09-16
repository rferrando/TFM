
import sys
import csv
from pprint import pprint
import numpy as np
import pandas as pd
from pandas import DataFrame
import itertools
import re

print("Reading Annotation Variants File")

maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)


file_target_genes = "/home/rosa/Master_Bioinformatica/TFM/scripts_Definitivos/Scripts/hgmd/Fármacos_priorizados_con_sus_dianas.csv"
sif_target_genes = "/home/rosa/Master_Bioinformatica/TFM/scripts_Definitivos/Scripts/hgmd/SIF_Fármacos_priorizados_con_sus_dianas.csv"
file_variants = "/home/rosa/Master_Bioinformatica/TFM/scripts_Definitivos/Scripts/hgmd/variants_to_add_to_nodes.csv"



def getting_genes_from_group(symbol):
    genes_set = set()
    if not symbol == None:
        symbol = symbol.strip('[]')
        list_genes = symbol.split("][")
#        print(list_genes)
        for gene in list_genes:
            genes_set.add(gene)
        return genes_set


def is_gene(s):
    pattern = '[^IES]'
    return re.match(pattern, s)


df_variants = pd.read_csv(file_variants, sep='\t')
df_sif = pd.read_csv(sif_target_genes, sep='\t')

genes = set()
i = 0

atts_variants = {}
sif_variants = {}
n = 0

with open(file_target_genes, "r") as t:
    reader = csv.reader(t,doublequote='False', quotechar='', quoting=csv.QUOTE_NONE)
    header = next(reader)  # skip the headers
    for node in reader:
        n = n + 1
        key = 'row' + str(n)
        other_node = {}
        other_node['name'] = node[3]
        other_node['Tipo Gen'] = node[0]
        other_node['symbol'] = node[1]
        other_node['symbol2'] = node[2]
        atts_variants.update({key: other_node})
        variant_node = {}
        if is_gene(node[0]):
            genes = set()
            genes.update(getting_genes_from_group(node[1]))
            df_variants_by_node = df_variants[df_variants['GENE'].isin(genes)]
            # if i==69: #'CYP3A5' in genes:
            #     print(i)
            #     print(genes)
            #     print(df_variants_by_node)
            for index, row in df_variants_by_node.iterrows():
                n = n + 1
                key = 'row' + str(n)
                variant_node = {}
                variant_node['name'] = row['HGMD']
                variant_node['Tipo Gen'] = 'V'
                if row['MUT'] == 'ALT':
                    variant_node['symbol'] = row['GENE'] + ': ' + row['ALT']
                    variant_node['symbol2'] = row['GENE'] + ': ' +row['ALT']
                else:
                    variant_node['symbol'] = row['GENE'] + ': ' +row['REF']
                    variant_node['symbol2'] = row['GENE'] + ': ' +row['REF']
                atts_variants.update({key: variant_node})
                # if i==69: #'CYP3A5' in genes:
                    # print(atts_variants)
                variant_node_sif = {}
                variant_node_sif['node1'] = variant_node['name']
                variant_node_sif['node2'] = other_node['name']
                variant_node_sif['interaction'] = 'variant interaction'
                sif_variants.update({key: variant_node_sif})
            i = i + 1

# print(i)
# print(len(genes))


# print(atts_variants)
df_network_nodes_variants =  pd.DataFrame.from_dict(atts_variants, orient='index', columns=header)
# df_priorizated_variants = df_variants[df_variants['GENE'].isin(genes)]
file_atts_variants_to_image = 'variants_for_image_atts.csv'
df_network_nodes_variants.to_csv(file_atts_variants_to_image, index=False, sep='\t')

df_network_interaction_variants =  pd.DataFrame.from_dict(sif_variants, orient='index', columns=['node1', 'interaction', 'node2'])

frames = [df_sif, df_network_interaction_variants]

result_df = pd.concat(frames)
file_sif_variants_to_image = 'variants_for_image_sif.csv'
result_df.to_csv(file_sif_variants_to_image, index=False, sep='\t')
