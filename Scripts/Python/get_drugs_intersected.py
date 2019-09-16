
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


file_atts = 'CriterioA-PG.csv'
Genes = []
header_atts = ""
with open(file_atts, "r") as a:
    reader = csv.reader(a)#delimiter='\t')
    # reader = csv.reader(a)
    header_atts = next(reader)  # skip the headers

    for row in reader:
        if row[12] == 'true':#'CriterioA-PG'
            attribute = []
            attribute = get_attribute("drugs", header_atts, row)
#			print(attribute)
            Genes.append({'Gene':attribute[0], 'drugs': attribute[1]})


Drugs = set()
def get_drugs(gene):
    for drug in gene['drugs'].split(","):
        if drug != "None":
            Drugs.add(drug)

for gene in Genes:
	get_drugs(gene)

file_drugs = 'Priorizated_drugs.csv'
df_drugs = pd.read_csv(file_drugs, sep='\t')
Drugs_priorizated = set()

for index, row in df_drugs.iterrows():
    Drugs_priorizated.add(row ['drug'])

intersected = Drugs.intersection(Drugs_priorizated)
# file = "Intersected_drugs.csv"
# file = open(file_sif, "w")

print("intersected drugs", len(intersected), intersected)
