
import csv
from pprint import pprint
import numpy as np
import pandas as pd
from pandas import DataFrame
import itertools

file_drugs = 'DrugHistogram_Info_output.csv'
df_drugs = pd.read_csv(file_drugs)#sep='\t'
df_drugs = df_drugs[df_drugs['number'] != 0]
df_drugs['list_of_targets'] = df_drugs.apply(lambda x: x['number_of_targets'][len(x['number_of_targets'].split(' ')[0]) + len(x['number_of_targets'].split(' ')[1]) + 2:], axis=1)


file_sif = "SIF_prueba_numero_primarios.csv"
df_sif = pd.read_csv(file_sif, sep='\t')



for index, row in df_drugs.iterrows():
    if row ['type'] == 'E':
        # print( row['drugs'], row['list_of_targets'].split(','))
        for gene in row['list_of_targets'].split(','):
            gene_info = {'Node 1': row['drugs'], 'Node 2': gene, 'interaction_type': 'drug_interaction', 'pathways': ""}
            # df = pd.DataFrame([{'c1':10, 'c2':100}]
            if not((df_sif['Node 1'] == row['drugs']) & (df_sif['Node 2'] == gene)).any():
                new_sif_df  = pd.DataFrame([gene_info], columns=df_sif.columns)
                df_sif = df_sif.append(new_sif_df, ignore_index=True)

file_sif_output = "SIF_prueba_numero_primarios_ampliado.csv"
df_sif.to_csv(file_sif_output)
