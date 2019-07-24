# -*- coding: utf-8 -*-
# MySQL Workbench Python script
# <description>
# Written in MySQL Workbench 6.0.8


#usr/share/doc
#"C:\Program Files\MySQL\MySQL Workbench 6.1 CE\MySQLWorkbench.exe" -query "Local instance MySQL56" -run-python "execfile('/home/practicas01/.mysql/workbench/scripts/hgmd_script.py')" -log-to-stderr -v

import grt
import csv
#import mforms
result = grt.root.wb.sqlEditors[0].executeScript("SELECT * from hgmd_pro.allmut where gene in (SELECT distinct(gene_sym) FROM hgmd_phenbase.hgmd_mutation JOIN hgmd_phenbase.hgmd_phenotype ON hgmd_mutation.phen_id = hgmd_phenotype.phen_id JOIN hgmd_phenbase.phenotype_concept ON hgmd_phenotype.phen_id = phenotype_concept.phen_id JOIN hgmd_phenbase.concept ON phenotype_concept.cui = concept.cui WHERE concept.str = 'endometriosis')")
resultset = result[0]
#print(type(resultset))

# Row Header
header = [i.name for i in resultset.columns]
#print(data)

# Data Rows
data = []
column_count = len(resultset.columns)
#print(column_count)
cursor = resultset.goToFirstRow()
#print(resultset.rowCount)

for n in range(1, resultset.rowCount):
    data.append([resultset.stringFieldValue(i) for i in range(column_count)])
    cursor = resultset.nextRow()
        
#print(data[64])
genes_DT = ""
with open("/home/practicas01/Desktop/rosa/hgmd/diseases.csv", "r") as csvfile:
    csv_reader = csv.reader(csvfile)
    genes_DT = next(csv_reader)
    genes_DT = "'" + str(genes_DT[0]).upper() + "'"
    for row in csv_reader:
        genes_DT = genes_DT + ',' + "'" + str(row[0]).upper() + "'"

with open("/home/practicas01/Desktop/rosa/hgmd/targets.csv", "r") as csvfile:
    csv_reader = csv.reader(csvfile)
    genes_DT_aux = next(csv_reader)
    genes_DT = genes_DT + "'" + str(genes_DT_aux[0]).upper() + "'"
    for row in csv_reader:
        genes_DT = genes_DT + ',' + "'" + str(row[0]).upper() + "'"

genes_clause = "SELECT * from hgmd_pro.allmut where gene in (" + genes_DT + ")"
#print(genes_clause)
result = grt.root.wb.sqlEditors[0].executeScript(genes_clause)

resultset = result[0]
#print(type(resultset))

# Data Rows
column_count = len(resultset.columns)
#print(column_count)
cursor = resultset.goToFirstRow()
#print(resultset.rowCount)

for n in range(1, resultset.rowCount):
    data.append([resultset.stringFieldValue(i) for i in range(column_count)])
    cursor = resultset.nextRow()

with open("/home/practicas01/Desktop/rosa/result.tsv", "w") as output_file:
    data_string = ""
    for item in header:
        data_string = data_string + item + '\t'
            
    data_string = data_string + '\n'
    output_file.write(data_string)
#    print(data_string)

    data_string = ""
    for index, item in enumerate(data):
        value_str = list(map(str, item))
        data_string = '\t'.join(value_str) + '\n'
        output_file.write(data_string)


