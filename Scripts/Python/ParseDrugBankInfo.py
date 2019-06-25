#!/usr/bin/env python

from utils.MapGeneSymbolToEntrez import Proxy
import csv
import numpy as np
import pandas as pd
from pandas import DataFrame

class DrugTargets:

    def __init__(self, file, columns, test_rows):
        self.df = pd.DataFrame(columns=columns)
        self.columns = columns
        self.proxy =  Proxy()
        self.test_rows = test_rows
        self.df_drugs = pd.read_csv(file, sep='\t')

        # df.dtypes # muestra los tipos de las columnas del DataFrame
        self._read_and_convert_file_to_drug_Target_info()

    def _is_nan(self, x):
        return (x is np.nan or x != x)

    def _split_genes(self, gene_text):
        return gene_text.split(";")

    def _get_all_target_genes_by_drug(self, drug_info, type_of_target):

        genes = []
        if not(self._is_nan(drug_info[type_of_target])):
            gs = self._split_genes(drug_info[type_of_target])
            gs = list(map(str.upper,gs))
            genes.extend(gs)

        #Remove duplicates
        genes = list(dict.fromkeys(genes))
        # print(genes)
        return genes

    def _copy_drug_info_and_add_ids(self, drug_info, gene_info, type_of_target):
        drug_target = drug_info[self.columns]

        drug_target['Gene'] = gene_info
        if type_of_target in ['transporters', 'enzymes', 'carriers']:
            drug_target['ADME'] = True

        return drug_target

    def _get_gene_info(self, drug_info, type_of_target):
        target_genes_of_drug = self._get_all_target_genes_by_drug(drug_info, type_of_target)
        # all_genes_ids = list(map(self.proxy.map_ids, target_genes_of_drug))

        for gene_ids in target_genes_of_drug:
            target_gene_info = self._copy_drug_info_and_add_ids(drug_info, gene_ids, type_of_target)
            target_gene_df  = pd.DataFrame([target_gene_info], columns=self.columns)
            self.df = self.df.append(target_gene_df , ignore_index=True)

        return target_genes_of_drug

    def _format_drug_target_info(self, drug_info):
        all_genes = []
# Obtenemos la info de los Targets y genes ADME que están aprobados
        targets_genes = self._get_gene_info(drug_info, "targets")
        transporters_genes = self._get_gene_info(drug_info, "transporters")
        enzymes_genes = self._get_gene_info(drug_info, "enzymes")
        carriers_genes = self._get_gene_info(drug_info, "carriers")

        #Join all kind of targets
        all_genes.extend(targets_genes)
        all_genes.extend(transporters_genes)
        all_genes.extend(enzymes_genes)
        all_genes.extend(carriers_genes)

        #Remove duplicates
        all_genes = list(dict.fromkeys(all_genes))
        #Get all of targets
        all_genes_str = ','.join(all_genes)

        return all_genes_str


    def _read_and_convert_file_to_drug_Target_info(self):
        if (self.test_rows) > 0:
            self.df_drugs = self.df_drugs[0:self.test_rows]

        self.df_drugs['all_targets'] = self.df_drugs.apply(self._format_drug_target_info, axis=1)

    def get_drug_info(self, file):
        self._write_df_to_file(file, self.df_drugs)
        return self.df_drugs

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
        df_group = self.df.groupby(['Gene'],as_index=False)['drugs','ID','ADME'].agg(lambda x: list(dict.fromkeys(x))) #(','.join)
        df_group.apply(self._update_vals, axis=1)
        DataFrame.drop_duplicates(df_group)
        self._write_df_to_file(file, df_group)
        return df_group

if __name__=="__main__":
    columns=['Gene', 'ADME', 'drugs', 'ID']
    print("Reading File with All Drugs")

    file = 'drugbank.csv'
    file_drugs_approved = 'drugbank_output.csv'
    output_file = open(file_drugs_approved, "w")
    with open(file, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        header = next(csv_reader)
        data = ""
        for value in header:
            data = data + str(value) + "\t"

        data = data + "\n"
        output_file.write(data)
        for row in csv_reader:
            list_groups = row[6].split(";")
            if "approved" in list_groups:
                data = ""
                for value in row:
                    data = data + str(value) + "\t"

                data = data + "\n"
                output_file.write(data)
    #
    print("Reading File with Approved Drugs")

    DrugTargets =  DrugTargets(file_drugs_approved, columns, 0)
    DrugInfo_df = DrugTargets.get_drug_info('DrugHistogram_Info.csv')
    DrugTargets_df = DrugTargets.get_drug_target_info('All_DrugBank_Info.csv')
    Drugs_by_gene_df = DrugTargets.get_drugs_by_gen('Drugs_by_gene_All_DruBank.csv')
