from utils.MapGeneSymbolToEntrez import Proxy
import csv
import numpy as np
import pandas as pd
from pandas import DataFrame
import functools
import itertools

# file_atts = "ATTS_directed_cluster_and_drugs_without_outliers_with_proximity_to_disease_genes.csv"
# file_atts = "ATTS_directed_cluster_and_drugs_only_with_proximity_to_disease_genes.csv"
# file_atts = "ATTS_directed_cluster_and_drugs_without_outliers_cluster22.csv"
# file_atts = "ATTS_directed_cluster_and_drugs.csv"
file_atts = "ATTS_SubNetwork_2Neighbors_filtered by diseases nodes OR nodes connected to diseases nodes(1path).csv"
df = pd.read_csv(file_atts)
# df = pd.read_csv(file_atts, sep='\t')
file_translate = 'symbol_identificators.csv'
df_translate = pd.read_csv(file_translate)#sep='\t'

disease_file = 'diseases.csv'
diseases = set()

with open(disease_file, "r") as d:
    csv_reader = csv.reader(d)
    for row in csv_reader:
        diseases.add(row[0])

target_file = 'targets.csv'
targets = set()

with open(target_file, "r") as t:
    csv_reader = csv.reader(t)
    for row in csv_reader:
        targets.add(row[0])

def get_all_genes_of_group(name):
    return name.split(' ')

def format_when_target_or_disease(symbol):
    if symbol in diseases or symbol in targets:
        symbol = '[' + symbol + ']'
    return symbol

def map_genes_of_group(total, gene):
#    print(total)
    and_group = []
    symbol = ""
    if "&" in gene:
        and_group = gene.split('&')
        try:
            symbol1 = df_translate[df_translate['entrezgene'] == int(and_group[0][4:len(and_group[0])])]['hgnc_symbol'].tolist()[0]
            symbol1 = format_when_target_or_disease(symbol1)
#        print(gene, symbol1)
        except IndexError as e:
            print(e, symbol1, and_group)
        try:
            symbol2 = df_translate[df_translate['entrezgene'] == int(and_group[1][4:len(and_group[1])])]['hgnc_symbol'].tolist()[0]
            symbol2 = format_when_target_or_disease(symbol2)
#        print(gene,symbol2)
        except IndexError as e:
            print(e, symbol2, and_group)

        symbol = symbol1 + '&' + symbol2
    else:
        try:
            symbol = df_translate[df_translate['entrezgene'] == int(gene[4:len(gene)])]['hgnc_symbol'].tolist()[0]
            symbol = format_when_target_or_disease(symbol)
#        print(gene, symbol)
        except IndexError as e:
            print(e, gene)

    return total + ' ' + symbol

def convert_node_to_symbol(atts):
    if atts['Tipo Gen'] in ['G','T','D','TD']:
        genes_of_group = get_all_genes_of_group(atts['name'])
#    print(genes_of_group)
        if len(genes_of_group) > 1 or ("&" in genes_of_group[0]):
            symbol = functools.reduce(map_genes_of_group, genes_of_group, '')
        else:
            gene = genes_of_group[0]
            try:
                symbol = df_translate[df_translate['entrezgene'] == int(gene[4:len(gene)])]['hgnc_symbol'].tolist()[0]
            except IndexError as e:
                print(e, gene)
    else:
        symbol = atts['symbol']
    return symbol


df['symbol2'] = df.apply(convert_node_to_symbol, axis=1)
# file_output = "ATTS_directed_cluster_and_drugs_without_outliers_with_proximity_to_disease_genes_with_all_node_names.csv"
# file_output = "ATTS_directed_cluster_and_drugs_only_with_proximity_to_disease_genes_with_all_node_names.csv"
# file_output = "ATTS_directed_cluster_and_drugs_without_outliers_cluster22_with_all_node_names.csv"
file_output = "ATTS_directed_cluster_and_drugs_with_all_node_names.csv"
file_output = "ATTS_SubNetwork_2Neighbors_filtered by diseases nodes OR nodes connected to diseases nodes(1path)_with_all_node_names.csv"
df.to_csv(file_output)
