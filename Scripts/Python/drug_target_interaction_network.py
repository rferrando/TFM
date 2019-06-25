
import csv
from pprint import pprint
import numpy as np
import pandas as pd
from pandas import DataFrame
import itertools

print("Reading File")

def get_attribute(attribute_name, header, row):
    gene_value = ""
    attribute_value = ""
    num_attr = len(row)
    indexes = list(range(num_attr))

    for i in indexes:
        if header[i] == "name":
            gene_value = str(row[i])
        if header[i] == attribute_name:
            attribute_value = str(row[i])

    return [gene_value, attribute_value]

file_drugs = 'DrugHistogram_Info_output.csv'
df = pd.read_csv(file_drugs)#sep='\t'
df = df[df['number'] != 0]

# file_atts = "ATTS_directed_cluster_and_drugs_without_outliers_with_proximity_to_disease_genes.csv"
# file_atts = "ATTS_directed_cluster_and_drugs_only_with_proximity_to_disease_genes.csv"
# file_atts = "ATTS_directed_cluster_and_drugs_without_outliers_cluster22.csv"
# file_atts = "ATTS_directed_cluster_and_drugs.csv"
file_atts = "ATTS_SubNetwork_2Neighbors_filtered by diseases nodes OR nodes connected to diseases nodes(1path).csv"
Genes = []
header_atts = ""
with open(file_atts, "r") as a:
    # reader = csv.reader(a, delimiter='\t')
    reader = csv.reader(a)
    header_atts = next(reader)  # skip the headers

    for row in reader:
        attribute = []
        attribute = get_attribute("drugs", header_atts, row)
#			print(attribute)
        Genes.append({'Gene':attribute[0], 'drugs': attribute[1]})

i = 0
Drugs = set()
def writing_SIF_for_node(gene, i):
    data = ""
    for drug in gene['drugs'].split(","):
        if drug != "None":
            # data =  data + drug + "\t" + "drug_interaction" + "\t" +  gene['Gene'] + "\t" + "" + "\n"
            data =  data + drug + "," + "drug_interaction"  + "," +  gene['Gene'] + "," + "" + ","+ "" + ","+ "" + "\n"
            Drugs.add(drug)
            i = i + 1
    return [data, i]

def writing_ATTS_for_node(header, node):
    data = ""
    num_attr = len(header)
    indexes = list(range(num_attr))

    for i in indexes:
        attribute_value = ""
        if header[i] == "name" or header[i] == "symbol":
            attribute_value = str(node)
        if header[i] == "Tipo Gen":
            try:
                attribute_value = df[df['drugs'] == str(node)]['type'].tolist()[0]
            except IndexError as e:
                print(e, str(node))
            # print(type(attribute_value))

        # data =  data + attribute_value + "\t"
        data =  data + attribute_value + ","

    data =  data + "\n"
    return data

# file_sif = "SIF_directed_cluster_and_drugs_without_outliers_with_proximity_to_disease_genes.csv"
# file_sif = "SIF_directed_cluster_and_drugs_only_with_proximity_to_disease_genes.csv"
# file_sif = "SIF_directed_cluster_and_drugs_without_outliers_cluster22.csv"
# file_sif = "SIF_directed_cluster_and_drugs.csv"
file_sif = "SIF_SubNetwork_2Neighbors_filtered by diseases nodes OR nodes connected to diseases nodes(1path).csv"
sif_file = open(file_sif, "a")
atts_file =open(file_atts, "a")

for gene in Genes:
	[data_SIF, i] = writing_SIF_for_node(gene, i)
	sif_file.write(data_SIF)

for drug in Drugs:
	data_ATTS = writing_ATTS_for_node(header_atts, drug)
	atts_file.write(data_ATTS)

print("interactions drug-target", i)
