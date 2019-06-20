#!/usr/bin/env python

from utils.MapGeneSymbolToEntrez import Proxy
import csv
import numpy as np
import pandas as pd
from pandas import DataFrame

class DrugTargets:

    def __init__(self, file, columns, test_rows):
        self.file = file
        self.df = pd.DataFrame(columns=columns)
        self.columns = columns
        self.proxy =  Proxy()
        self.test_rows = test_rows
        self.drugs_T = set()

        self._read_and_convert_file_to_drug_Target_info()

    def _is_nan(self, x):
        return (x is np.nan or x != x)

    def _split_genes(self, gene_text):
        return gene_text.split(";")

    def _get_all_target_genes_by_drug(self, drug_info, type_of_target):

        genes = []
        try:
            if not(self._is_nan(drug_info[type_of_target])):
                gs = self._split_genes(drug_info[type_of_target])
                gs = list(map(str.upper,gs))
                genes.extend(gs)

            #Remove duplicates
            genes = list(dict.fromkeys(genes))
            # print(genes)
        except KeyError as e:
            print(drug_info['drugs'], e)
        return genes

    def _copy_drug_info_and_add_ids(self, drug_info, gene_info, type_of_target):
        drug_target = drug_info[self.columns]

        try:
            drug_target['Gene ID'] = gene_info['entrezgene']
        except KeyError:
            drug_target['Gene ID'] = 'NaN'
        drug_target['Gene'] = gene_info['symbol']

        if type_of_target in ['transporters', 'enzymes', 'carriers']:
            drug_target['ADME'] = True

        return drug_target

    def _get_gene_info(self, drug_info, type_of_target):
        target_genes_of_drug = self._get_all_target_genes_by_drug(drug_info, type_of_target)
        all_genes_ids = list(map(self.proxy.map_ids, target_genes_of_drug))

        for gene_ids in all_genes_ids:
            target_gene_info = self._copy_drug_info_and_add_ids(drug_info, gene_ids, type_of_target)
            # if target_gene_info['Gene'] =='BCL10':
            # if drug_info['drugs'] == 'Benzoic Acid' or drug_info['drugs'] ==' B' or drug_info['drugs'] ==' A' or drug_info['drugs'] =='Ganciclovir' or drug_info['drugs'] =='Synthetic Conjugated Estrogens A' or drug_info['drugs'] =='Synthetic Conjugated Estrogens' or drug_info['drugs'] =='Interferon Alfa-2b Recombinant' or drug_info['drugs'] =='Synthetic Conjugated Estrogens B':
            #
            #     print('drug: ', drug_info['ID'])
            #     print('info: ', target_genes_of_drug)
            #     print('tipo: ',type_of_target)
            #     print('Genes: ',all_genes_ids)
            target_gene_df  = pd.DataFrame([target_gene_info], columns=self.columns)
            self.df = self.df.append(target_gene_df , ignore_index=True)

    def _format_drug_target_info(self, drug_info):
# Obtenemos los IDs de los Targets
        self._get_gene_info(drug_info, "targets")
        self._get_gene_info(drug_info, "transporters")
        self._get_gene_info(drug_info, "enzymes")
        self._get_gene_info(drug_info, "carriers")

        if drug_info['ID'] != None:
            self.drugs_T.update(drug_info['ID'].split(','))

    def _read_and_convert_file_to_drug_Target_info(self):
        df = pd.read_csv(self.file, sep='\t', dtype={'transporters':'str','targets':'str','enzymes':'str','carriers': 'str'})
        print(df.dtypes)# muestra los tipos de las columnas del DataFrame
        if (self.test_rows) > 0:
            df = df[0:self.test_rows]

        df.apply(self._format_drug_target_info, axis=1)
        file_drugs_T = open('drugs_T.csv', 'w') # open for 'w'riting
        file_drugs_T.write("\n".join(self.drugs_T))

    def get_drug_target_info(self, file):
        DataFrame.drop_duplicates(self.df)
        self._write_df_to_file(file, self.df)
        return self.df

    def _write_df_to_file(self, file, df):
        df.to_csv(file)

    def _clean_data(self, row, field):

        try:
            data = row[field]
        except AttributeError as e:
            pass

        for i,item in enumerate(data):
            data[i]=str(item)
        try:
            data.remove('nan')
        except ValueError:
            pass
        data = list(dict.fromkeys(data))
        data = ','.join(data)
        # print(row[field], data)
        return data

    def _update_vals(self, row):
        row['drugs'] = self._clean_data(row, 'drugs')
        row['ID'] = self._clean_data(row, 'ID')
        row['ADME'] = self._clean_data(row, 'ADME')

        return row


    def get_drugs_by_gen(self, file):
        df_group = self.df.groupby(['Gene', 'Gene ID'],as_index=False)['drugs','ID','ADME'].agg(lambda x: list(dict.fromkeys(x))) #(','.join)
        df_group.apply(self._update_vals, axis=1)
        DataFrame.drop_duplicates(df_group)
        self._write_df_to_file(file, df_group)
        return df_group

if __name__=="__main__":
    columns=['Gene', 'Gene ID', 'ADME', 'drugs', 'ID']

    file = 'approved_drugs_endometriosis_reduced.csv'
    DrugTargets =  DrugTargets(file, columns, 0)
    DrugTargets_df = DrugTargets.get_drug_target_info('Drug_Targets_Info_TA.csv')
    Drugs_by_gene_df = DrugTargets.get_drugs_by_gen('Drugs_by_gene_TA.csv')
