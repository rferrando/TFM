#!/usr/bin/env python

from utils.MapGeneSymbolToEntrez import Proxy
import csv
import numpy as np
import pandas as pd
from pandas import DataFrame

class GWAS:

    def __init__(self, file, columns, test_rows):
        self.file = file
        self.GWAS_df = pd.DataFrame(columns=columns)
        self.columns = columns
        self.proxy =  Proxy()
        self.test_rows = test_rows

        self._read_and_convert_file_to_GWAS_variant_info()


    def _split_genes(self, gene_text):
        _split_genes = []
        for text in gene_text.split(", "):
            _split_genes.extend(text.split(" - "))

        return _split_genes

    def _is_nan(self, x):
        return (x is np.nan or x != x)

    def _get_all_genes_for_variant_by_study(self, GWAS_variant_info):

        genes = []
        if (GWAS_variant_info['Gene'] != 'NR'
        and GWAS_variant_info['Gene'] != 'intergenic'
        and GWAS_variant_info['Gene'] != 'Intergenic'
        and not(self._is_nan(GWAS_variant_info['Gene']))):
            gs = self._split_genes(GWAS_variant_info['Gene'])
            genes.extend(gs)

        if not(self._is_nan(GWAS_variant_info['Gene 2'])):
            gs = self._split_genes(GWAS_variant_info['Gene 2'])
            genes.extend(gs)

        if not(self._is_nan(GWAS_variant_info['OrherGenes'])):
            gs = self._split_genes(GWAS_variant_info['OrherGenes'])
            genes.extend(gs)

        #Remove duplicates
        genes = list(dict.fromkeys(genes))
        return genes

    def _copy_variant_info_and_add_ids(self, GWAS_variant_info, gene_info):
        study_variant = GWAS_variant_info[self.columns]

        try:
            study_variant['Gene ID'] = gene_info['entrezgene']
        except KeyError:
            study_variant['Gene ID'] = 'NaN'
        study_variant['Gene'] = gene_info['symbol']
        return study_variant

    def _format_GWAS_info(self, GWAS_variant_info):
        genes_of_variant = self._get_all_genes_for_variant_by_study(GWAS_variant_info)
        all_genes_ids = list(map(self.proxy.map_ids, genes_of_variant))

        for gene_ids in all_genes_ids:
            gene_variant_info = self._copy_variant_info_and_add_ids(GWAS_variant_info, gene_ids)
            gene_variant_df  = pd.DataFrame([gene_variant_info], columns=self.columns)
            self.GWAS_df = self.GWAS_df.append(gene_variant_df , ignore_index=True)

    def _read_and_convert_file_to_GWAS_variant_info(self):
        GWAS_df = pd.read_csv(self.file, dtype={'SNP rs':'str','Gene ID':'str','Gene ID 2':'str','Chromosome': 'str'})
        # GWAS_df.dtypes # muestra los tipos de las columnas del DataFrame
        if (self.test_rows) > 0:
            GWAS_df = GWAS_df[0:self.test_rows]

        GWAS_df.apply(self._format_GWAS_info, axis=1)

    def get_whole_GWAS(self, file):
        self._write_GWAS_to_file(file, self.GWAS_df)
        return DataFrame.drop_duplicates(self.GWAS_df)

    def _write_GWAS_to_file(self, file, df):
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
        fields = ['Gene', 'Gene ID', 'Population','PubMed','Context', 'P-Value']
        for i,field in enumerate(fields):
            row[field] = self._clean_data(row, field)
        return row

    def _update_vals2(self, row):
        fields = ['SNP rs','Chromosome', 'Location', 'Population','PubMed','Context', 'P-Value']
        for i,field in enumerate(fields):
            row[field] = self._clean_data(row, field)
        return row

    def get_genes_by_SNP(self, file):
        df_group = self.GWAS_df.groupby(['SNP rs','Chromosome', 'Location'],as_index=False)['Gene', 'Gene ID', 'Population','PubMed','Context', 'P-Value'].agg(lambda x: list(dict.fromkeys(x))) #(','.join)
        df_group.apply(self._update_vals, axis=1)
        DataFrame.drop_duplicates(df_group)
        self._write_GWAS_to_file(file, df_group)
        return df_group
    def get_SNPs_by_gen(self, file):
        df_group = self.GWAS_df.groupby(['Gene', 'Gene ID'],as_index=False)['SNP rs','Chromosome', 'Location', 'Population','PubMed','Context', 'P-Value'].agg(lambda x: list(dict.fromkeys(x))) #(','.join)
        df_group.apply(self._update_vals2, axis=1)
        DataFrame.drop_duplicates(df_group)
        self._write_GWAS_to_file(file, df_group)
        return df_group

if __name__=="__main__":
    GWAS_file = 'GWAS_Info.csv'
    GWAS_columns=['PubMed', 'SNP rs', 'Gene', 'Gene ID', 'Chromosome', 'Location', 'Context', 'P-Value', 'Source', 'Population']
    #GWAS_columns=['PubMed', 'Trait', 'SNP rs', 'Context', 'Gene', 'Gene ID', 'Chromosome', 'Location', 'P-Value', 'Source', 'Population']
    GWAS =  GWAS(GWAS_file, GWAS_columns, 0)

    Whole_GWAS_df = GWAS.get_whole_GWAS('Whole_GWAS_formatted.csv')
    # print(Whole_GWAS_df)

    GWAS_by_gene_df = GWAS.get_SNPs_by_gen('GWAS_by_gene_formatted.csv')
    # print(GWAS_by_gene_df)

    GWAS_by_SNP_df = GWAS.get_genes_by_SNP('GWAS_by_SNP_formatted.csv')
    # print(GWAS_by_SNP_df)
