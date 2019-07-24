
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

file_variants = 'variants_with_high_impact_regardless_rankscore.csv'
file_target_genes = 'genes_dianas_de_drogas_con_outdegree>30_and_8.csv'
df_variants = pd.read_csv(file_variants, sep='\t')

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

genes = set()
i = 0
with open(file_target_genes, "r") as t:
    reader = csv.reader(t, delimiter='\t',doublequote='False', quotechar='', quoting=csv.QUOTE_NONE)
    header = next(reader)  # skip the headers
    # header = ['ADME,Analisis_Variantes,cl1.Status,clinvar_clnsig,Cluster Number,Context,dbsnp,disease,drugs,mutype,name,Nº Variantes,NºDrogas,primary,symbol,symbol2,tag,target,Tipo Gen']
    for row in reader:
        if is_gene(row[18]): #and row[0] and row[2] != 'Outlier':
            genes.update(getting_genes_from_group(row[14]))
            i = i + 1

print(i)
print(genes)
#df_variants = df_variants[df_variants['RANKSCORE'] > 0.98]
#['RANKSCORE'] > 0.98 (for all genes- 51)--> (208, 18)
#['RANKSCORE'] > 0.98 (for ADME genes- 6)--> (74, 18)
#['RANKSCORE'] > 0.98 (for ADME genes in clusters- 1)--> (4, 18)

df_priorizated_variants = df_variants[df_variants['GENE'].isin(genes)]
#df_priorizated_variants = df_priorizated_variants[df_priorizated_variants['CLASS'].isin(['DM', 'DM?', 'DP', 'DFP'])]
# df_priorizated_variants = df_priorizated_variants[df_priorizated_variants['PHEN'].str.contains('endometriosis')]
type_variant_set = set()
def get_type_variant(transcript):
    transcript_list = transcript.split('|')
    if transcript_list[2] == 'HIGH':
        return transcript_list[1]
    else:
        return ""

def parse_ANN_info(info):
#Obtener el impacto de todos los transcritos
    ANN = info['ANN'].split(',')
    map_iterator = map(get_type_variant, ANN)
    ANN = list(map_iterator)
    type_variant_set .update(ANN)
    # attributes['IMPACT'] = ['MODIFIER', 'MODIFIER', 'MODIFIER', 'MODIFIER', 'MODIFIER', 'MODIFIER', 'MODIFIER']
    return dict(zip(ANN,map(ANN.count,ANN)))

df_priorizated_variants['variant_type'] = df_priorizated_variants.apply(parse_ANN_info, axis=1)

print(df_priorizated_variants.iloc[0])
print(df_priorizated_variants.shape)
print("\n".join(type_variant_set))

file_variants_to_analyse = 'variants_to_analyse(outdegree>30_and_8).csv'
df_priorizated_variants.to_csv(file_variants_to_analyse, index=False, sep='\t')
