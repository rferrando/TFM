#Integra la información genética de los genes Targets y Diseases en la subred generada
from igraph import *
import csv
import argparse
import re
import itertools
from pprint import pprint

parser = argparse.ArgumentParser(description='Create SubNetwork from interest paths')
parser.add_argument("--interactome", help="KEGG SIF file", required = True)
parser.add_argument("--SIF", help="Output SIF file", required = True)
parser.add_argument("--ATTS", help="KEGG ATTS file", required = True)
parser.add_argument("--neighbors", help="Number of neighbors", type=int, default = 1)
parser.add_argument("--target-genes", help="Target Genes", required = True)
parser.add_argument("--disease-genes", help="Disease genes", required = True)

args = parser.parse_args()

#KEGG Network
g = Graph(directed=True)
nodes = set()
arcs = []
#Edge attributes
interaction_type = []
pathways = []

#Target and Disease
targets = []
diseases = []
sif_header = []

print("Creating graph")

with open(args.interactome, 'r') as csvfile:
	csv_reader = csv.reader(csvfile)
	next(csv_reader)
	for row in csv_reader:
#Nodes
		nodes.add(row[0])
		nodes.add(row[1])
#Edges
		arcs.append((row[0], row[1]))
#Edge Attributes
		interaction_type.append(row[2])
		pathways.append(row[3])


for node in nodes:
	g.add_vertices(node)

g.add_edges(arcs)
g.es["interaction_type"] = interaction_type
g.es["pathways"] = pathways

def remove_duplicates_from_a_list(mylist):
  return list(dict.fromkeys(mylist))

def add_attributes(attributes, header, row, from_column):
	attr = {}
	num_attr = len(row)
	indexes = list(range(num_attr))
	for _ in itertools.repeat(None, from_column):
		indexes.pop(0) #Remove first columns
	for i in indexes:
		try:
			attr[header[i]] = attributes[header[i]] + ',' + row[i]
			attr_list = attr[header[i]].split(",")
			attr_list = remove_duplicates_from_a_list(attr_list)
			attr[header[i]] = ','.join(attr_list)
		except Exception as e:
			# print(e)
			attr[header[i]] = row[i]

	return attr
# print("Reading Genes ADME")
#
# ADME_file = 'genesAdme.txt'
# ADME_genes = set()
# with open(ADME_file, "r") as ADME:
# 	csv_reader = csv.reader(ADME)
# 	next(csv_reader)
# 	for row in csv_reader:
# 		ADME_genes.add(row[0].upper())

print("Reading Drugs of targets")
drugs_of_targets_file = 'drugs_T.csv'
drugs_T = set()

with open(drugs_of_targets_file, "r") as d:
    csv_reader = csv.reader(d)
    for row in csv_reader:
        drugs_T.add(row[0])

print("Reading HGMD Variants")

HGMD_file = 'Variants_by_gene_formatted.csv'
HGMD_genes = []

with open(HGMD_file, "r") as HGMD:
	csv_reader = csv.reader(HGMD)
	HGMD_header = next(csv_reader)
	for row in csv_reader:
		variants_attr = {}
		variants_attr = add_attributes(variants_attr, HGMD_header, row, 1)
		# if row[0] == 'CFTR':
		# 	print(variants_attr)
		HGMD_genes.append({'Gene':row[0], 'Atts': variants_attr})

#print(HGMD_genes)

def get_gene_info(list , mygen):
    return [d['Atts'] for d in list if d['Gene'] == mygen]

print("Reading DrugBank Drugs")

DrugBank_file = 'Drugs_by_gene_All_DruBank.csv'
DrugBank_genes = []

with open(DrugBank_file, "r") as DrugBank:
	csv_reader = csv.reader(DrugBank)
	DrugBank_header = next(csv_reader)
	for row in csv_reader:
		drugs_attr = {}
		drugs_attr = add_attributes(drugs_attr, DrugBank_header, row, 1)
		DrugBank_genes.append({'Gene':row[0], 'Atts': drugs_attr})

print("Reading Primary Targets")
endometriosis_drugs = set()
primary_targets_file = 'seed_drugs_targets_endometriosis.csv'
primary_targets = set()
with open(primary_targets_file, "r") as primary:
	csv_reader = csv.reader(primary)
	next(csv_reader)
	for row in csv_reader:
		endometriosis_drugs.add(row[0])
		primary_targets.update(list(map(str.strip, str(row[3]).split(';')))) #Targets
		primary_targets.update(list(map(str.strip, str(row[2]).split(';')))) #Transporters
		primary_targets.update(list(map(str.strip, str(row[4]).split(';')))) #Enzymes
		primary_targets.update(list(map(str.strip, str(row[5]).split(';')))) #Carriers

primary_targets = list(filter(None, primary_targets)) # remove empty elements
no_primary_targets = len(primary_targets)

print("Reading Targets and Diseases and all its attributes")
drugs_E = set()
def parse_file(file, type):
	DT_genes_mapped_KEGG = set()
	number_of_genes = 0
	no_mapped_kegg = 0
	with open(file, "r") as attributes_file:
		csv_reader = csv.reader(attributes_file)
		header = next(csv_reader)
		for row in csv_reader:
			attributes_parsed = {}
			if not(row[1] == 'NaN'):
				number_of_genes = number_of_genes + 1
				kegg_id = 'hsa:'+ row[1]
				pattern = re.compile(r'\b{}\b'.format(kegg_id))
				kegg_ids_within_groups = []
				kegg_ids_within_groups = g.vs.select(lambda vertex: bool(pattern.search(vertex['name'])))

				if (len(kegg_ids_within_groups) > 0):
					DT_genes_mapped_KEGG.add(row[0])
					for vertex in kegg_ids_within_groups:
						group= len(vertex['name'].split("hsa:")) - 1
						if (group > 1):
							try:
								g.vs[vertex.index]['symbol'] = g.vs[vertex.index]['symbol'] + '[' + row[0] + ']'
							except Exception as e:
								# print(e)
								g.vs[vertex.index]['symbol'] = '[' + row[0] + ']'
						else:
							g.vs[vertex.index]['symbol'] = row[0]
						attributes_parsed = add_attributes(g.vs[vertex.index], header, row, 2)
						for key, value in attributes_parsed.items():
							g.vs[vertex.index][key] = value
						# if vertex['name'] == 'hsa:100861540 hsa:1551 hsa:1577':
						# 	print(type, kegg_id, row[0], g.vs[vertex.index]['drugs'])
						# if (row[0] in ADME_genes):
						# 	g.vs[vertex.index]['ADME'] = True
						if "T" in type:
							g.vs[vertex.index]['target'] = True
							# Informo si el Target es primario o no
							if (row[0] in primary_targets):
								g.vs[vertex.index]['primary'] = True
						if "D" in type:
							g.vs[vertex.index]['disease'] = True
							# Informo las drogas que lo tienen como TARGET
							drugs = get_gene_info(DrugBank_genes, row[0])
							if drugs:
								drugs_E.update(drugs[0]['ID'].split(','))
								attributes_parsed  = add_attributes(g.vs[vertex.index], list(drugs[0].keys()), list(drugs[0].values()), 0)
								for key, value in attributes_parsed.items():
									g.vs[vertex.index][key] = value
								if vertex['name'] == 'hsa:1436 hsa:1956 hsa:1969 hsa:2260 hsa:2261 hsa:2263 hsa:2264 hsa:2321 hsa:2324 hsa:3480 hsa:3643 hsa:3791 hsa:3815 hsa:4233 hsa:4804 hsa:5156 hsa:5159 hsa:7010':
									print(type, kegg_id, row[0], g.vs[vertex.index]['drugs'])
						#Informo las variantes de HGMD
						variants = get_gene_info(HGMD_genes, row[0])
						if variants:
				#			if row[0] == 'ADRA1A' or row[0] == 'ADRA1B' or row[0] == 'ADRA1D' or row[0] == 'ADRA2A' or row[0] == 'ADRA2B':
							attributes_parsed  = add_attributes(g.vs[vertex.index], list(variants[0].keys()), list(variants[0].values()), 0)
							for key, value in attributes_parsed.items():
								g.vs[vertex.index][key] = value
				else:
					no_mapped_kegg = no_mapped_kegg + 1

	return [DT_genes_mapped_KEGG, (number_of_genes - no_mapped_kegg), number_of_genes]

nodes_u_targets = parse_file(args.target_genes, "T")
nodes_u_diseases = parse_file(args.disease_genes, "D")

print("Calculating paths")

def find_all_paths(graph, start, path_length, correspondence_id_node):
    def find_all_paths_aux(adjlist, start, path_length, path, correspondence_id_node):
        path = path + [start.index]

        paths = []
        if len(path) == path_length:
        	if start.attributes()["target"] or start.attributes()["disease"]:
        		return [path]
        	else:
        		return []

        elif len(path) >= 2:
        	if start.attributes()["target"] or start.attributes()["disease"]:
        		return [path] # Revisar cuando devuelve paths de tamano menor del total.....

        for node in adjlist[start.index] - set(path):
        	neighbor_node = correspondence_id_node[node]

        	if neighbor_node.attributes()["target"] or neighbor_node.attributes()["disease"]:
        		# print("-" + str([path + [node]]))
        		paths.extend([path + [node]])
        	else:
        		paths.extend(find_all_paths_aux(adjlist, neighbor_node, path_length, path, correspondence_id_node))
        return paths

    adjlist = [set(graph.neighbors(node)) for node in range(graph.vcount())]
    return find_all_paths_aux(adjlist, start, path_length, [], correspondence_id_node)

paths = []
correspondence = dict()
correspondence_name_id = dict()
correspondence_id_node = dict()

for node in g.vs:
	correspondence[str(node.index)] = node.attributes()["name"]
	correspondence_name_id[node.attributes()["name"]] = node
	correspondence_id_node[node.index] = node

# pprint(correspondence)
for node in g.vs:
	if not node.attributes()["target"] and not node.attributes()["disease"]:
		continue

	gene = node.attributes()['name']
	gene_paths = find_all_paths(g, correspondence_name_id[gene], path_length = args.neighbors + 2, correspondence_id_node = correspondence_id_node)

	for gene_path in gene_paths:
		paths.extend([gene_path])

def check_target_is_secondary(node):
	if node.attributes()["disease"]:
		return False #Si el nodo es target y disease a la vez
	else:
		if node.attributes()["target"]:
			if not(node.attributes()["primary"]):
				return True #Es Target secundario
			else: False #Es Target primario
		else: True #Es un intermediario que a efectos de no relevancia es como un target secundario


def check_path_connects_two_secondary_targets(new_path):
	first_node_in_path = g.vs.find(name=new_path[0])
	last_node_in_path = g.vs.find(name=new_path[len(new_path) - 1])
	if check_target_is_secondary(first_node_in_path):
		if check_target_is_secondary(last_node_in_path):
			return True
		else: False
	else: False

#Convert paths to gene names & filter interest paths
paths_hsa = list()

for p in paths:
	new_path = list()

	for elem in p:
		new_path.append(correspondence_id_node[elem].attributes()["name"])
	# paths_hsa.extend([new_path])
#Check if Path connects 2 Secondary Targets
	if not(check_path_connects_two_secondary_targets(new_path)):
		paths_hsa.extend([new_path])

print("Removing duplicates and Writing SIF y ATTS")
output_file = open(args.SIF, "w")
atts_file = open(args.ATTS, "w")

target_drugs = set()
disease_drugs = set()

def get_edge_tuple_with_ids(arc):
	return (g.vs[arc.tuple[0]]['name'],g.vs[arc.tuple[1]]['name'])

def writing_SIF_for_edge(arc):
	data = ""
	header = ""

	edge_tuple = get_edge_tuple_with_ids(arc)
	sif_attributes = arc.attributes()
	data = edge_tuple[0] + "\t"+ edge_tuple[1] + "\t"

	header = 'Node 1' + "\t"+ 'Node 2' + "\t"
	for key, value  in sif_attributes.items():
		data = data + str(value) + "\t"
		header = header + key + "\t"

	data = data + "\n"
	header = header + "\n"
	return [data, header]

def adding_node_attributes(node, node_attributes):
	data = ""
	header = ""
	num_dbsnp = 0
#Añadimos nº de drogas
	try:
		drugs_list = node_attributes['drugs'].split(",")
		drugs_list = remove_duplicates_from_a_list(drugs_list)
		num_drugs = len(drugs_list)
	except AttributeError as e:
		num_drugs = 0

	data = data + str(num_drugs) + "\t"
	header = header + "NºDrogas" + "\t"
#Añadimos nº variantes
	try:
		dbsnp_list = node_attributes['dbsnp'].split(",")
		dbsnp_list = remove_duplicates_from_a_list(dbsnp_list)
		num_dbsnp = len(dbsnp_list)
	except AttributeError as e:
		num_dbsnp = 0

	data = data + str(num_dbsnp) + "\t"
	header = header + "Nº Variantes" + "\t"
#Añadimos Tipo de Gen
	type = "G"
	if node_attributes['target']:
		target_drugs.add(node_attributes['name'])
		type = "T"
	if node_attributes['disease']:
		if type =="T":
			type = type + "D"
		else:
			type = "D"
			disease_drugs.add(node_attributes['name'])

	data = data + str(type) + "\t"
	header = header + "Tipo Gen" + "\t"

	data = data + "\n"
	header = header + "\n"
	return [data, header]


def writing_ATTS_for_node(node):
	data = ""
	header = ""
	node_attributes = [n for n in node][0].attributes()

	for key, value in node_attributes.items():
		data = data + str(value) + "\t"
		header = header + key + "\t"
	[data2, header2] = adding_node_attributes(node, node_attributes)
	data = data + data2
	header = header + header2
	return [data, header]

def getting_gene_symbol_for_intersection(node):
    symbol = node[0]['symbol']
    if not symbol == None:
        symbol = symbol.strip('[]')
        list_genes = symbol.split("][")
#        print(list_genes)
        for gene in list_genes:
            atts.add(gene)



sif = set()
atts = set()
sx_network = set()
header_time = True

for path in paths_hsa:
	for e in range(0, len(path) - 1):
# Removing duplicates
		if path[e] + "\t" + path[e + 1] in sif or path[e + 1] + "\t" + path[e] in sif:
			continue
		sif.add(path[e] + "\t" + path[e + 1])

#Getting edge, source_node and target_node from every path
		source_node = g.vs.select(name=path[e])
		target_node = g.vs.select(name=path[e + 1])
		edge = g.es.select(_between=(source_node, target_node))

# Getting gene symbol in subnetwork for intersection (with targets and diseases genes)
		getting_gene_symbol_for_intersection(source_node)
		getting_gene_symbol_for_intersection(target_node)
# Writing SIF for edge
		i = 0
		for arc in edge:
			[data_file, header_file] = writing_SIF_for_edge(arc)
#Writing SIF header
			if header_time and i == 0:
				output_file.write(header_file)
			i = i + 1
#Writing Edge Information
			output_file.write(data_file)
# Writing ATTS for source_node
		if not path[e] in sx_network:
			[data_file, header_file] = writing_ATTS_for_node(source_node)
			sx_network.add(path[e])
#Writing ATTS header
			if header_time:
				atts_file.write(header_file)
			atts_file.write(data_file)
			header_time = False
# Writing ATTS for target_node
		if not path[e + 1] in sx_network:
			[data_file, header_file] = writing_ATTS_for_node(target_node)
			sx_network.add(path[e + 1])
			atts_file.write(data_file)

print("Analyzing Network")
print("Targets primarios de drogas de Endometriosis", no_primary_targets)
print("Targets intersection with KEGG: ", nodes_u_targets[1], ' (primary ', len(nodes_u_targets[0].intersection(primary_targets)),')', ' de ', nodes_u_targets[2])
print("Diseases intersection with KEGG: ", nodes_u_diseases[1], ' de ', nodes_u_diseases[2])
print("KEGG Targets in SubNetwork: ", len(atts.intersection(nodes_u_targets[0])), ' (primary ', len(atts.intersection(primary_targets)),')', ' de ', nodes_u_targets[1])
file_targets = open('targets.csv', 'w') # open for 'w'riting
file_targets.write("\n".join(nodes_u_targets[0]))
print("KEGG Diseases in SubNetwork: ", len(atts.intersection(nodes_u_diseases[0])), ' de ', nodes_u_diseases[1])
print(atts.intersection(nodes_u_diseases[0]))
file_disease = open('diseases.csv', 'w') # open for 'w'riting
file_disease.write("\n".join(nodes_u_diseases[0]))
drugs_E_file = open('drugs_E.csv', 'w') # open for 'w'riting
drugs_E_file.write("\n".join(drugs_E))
print("KEGG Diseases not in SubNetwork: ", nodes_u_diseases[0].difference(atts))


print("Drug Analysis")
import numpy as np
import pandas as pd
from pandas import DataFrame
##############################
file_translate = 'symbol_identificators.csv'
df_translate = pd.read_csv(file_translate)#sep='\t'
##############################

def gene_name_is_in(drug_targets, genes_set, type, genes):
	i = 0
	for gene in drug_targets:
		try:
			entrezgene = df_translate[df_translate['hgnc_symbol'] == gene]['entrezgene'].tolist()[0]
			kegg_id = 'hsa:'+ str(entrezgene)
		except IndexError as e:
			# print(e, gene)
			continue

		for kegg_gene in genes_set:
			pattern = re.compile(r'\b{}\b'.format(kegg_id))
			if bool(pattern.search(kegg_gene)):
				i = i + 1
				genes.add(kegg_gene)

	return [i, genes]


def get_number_of_targets(drug_targets, set1=endometriosis_drugs, set2=target_drugs, set3=disease_drugs):
	genes = set()
	origin = set()
	type = ''

	gene_name_is_in(str(drug_targets['all_targets']).split(','), set2, 'T', genes)
	gene_name_is_in(str(drug_targets['all_targets']).split(','), set3, 'D', genes)

	if drug_targets['ID'] in endometriosis_drugs:
		type = 'S'
	elif drug_targets['ID'] in drugs_T:
	    	type = 'I'
	elif drug_targets['ID'] in drugs_E:
		type = 'E'

	# if drug_targets['drugs'] == 'Trastuzumab emtansine':
	# 	print (genes)
	# 	print (len(genes))
	return str(len(genes)) + ' ' + type + ' ' + ','.join(genes)

file_drugs = 'DrugHistogram_Info.csv'
df = pd.read_csv(file_drugs)#sep='\t'
df['number_of_targets'] = df.apply(get_number_of_targets, axis=1)
df['number'] = df.apply(lambda x:float(x['number_of_targets'].split(' ')[0]), axis=1)
df['type'] = df.apply(lambda x: x['number_of_targets'].split(' ')[1], axis=1)

file_drugs_output = 'DrugHistogram_Info_output.csv'
df.to_csv(file_drugs_output)

print("Drug Boxplot")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

df = df[df['number'] != 0]
data_plot = [df[df['type'] == 'E']['number'], df[df['type'] == 'I']['number'], df[df['type'] == 'S']['number']]
fig, ax = plt.subplots(figsize=(20, 10))
bp = ax.boxplot(data_plot)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')
colors = ['darkkhaki', 'royalblue', 'salmon']

ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

# Hide these grid behind plot objects
ax.set_axisbelow(True)
ax.set_title('Comparación del nº de targets para cada tipo de droga')
ax.set_xlabel('Tipo de droga')
ax.set_ylabel('Nº de targets')
num_boxes = len(data_plot)
for i in range(num_boxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    box_coords = np.column_stack([boxX, boxY])
    # Alternate between Dark Khaki and Royal Blue
    ax.add_patch(Polygon(box_coords, facecolor=colors[i % 3]))
# Multiple box plots on one Axes
ax.set_xticklabels(['Drogas cuyos targets son genes causantes de endometriosis','Drogas que interactuan con drogas de endometriosis', 'Drogas relacionadas con endometriosis'],
                    fontsize=9)#rotation=45, fontsize=8)
# Calculate number of obs per group & median to position labels
medians = df.groupby(['type'])['number'].median().sort_index(axis=0).values
nobs = df['type'].value_counts().sort_index(axis=0).values
nobs = [str(x) for x in nobs.tolist()]
nobs = ["n: " + i for i in nobs]
# nobs.reverse()
# Add it to the plot
pos = range(len(nobs))
for tick,label in zip(pos,ax.get_xticklabels()):
    ax.text(pos[tick] + 1, medians[tick] + 0.03, nobs[tick],
    horizontalalignment='center', size='x-small', color='w', weight='semibold')

plt.savefig('Drug_Boxplot.svg')
#plt.show()

print("Drug Density Plot")
import seaborn as sns
# plot of 2 variables
plt.figure()
p1=sns.kdeplot(df[df['type'] == 'S']['number'], shade=True, color=colors[1], label ='Drogas relacionadas con endometriosis')
p2=sns.kdeplot(df[df['type'] == 'I']['number'], shade=True, color=colors[2],label ='Drogas que interactuan con drogas de endometriosis')
p3=sns.distplot(df[df['type'] == 'E']['number'], color=colors[0], label ='Drogas cuyos targets son genes causantes de endometriosis', hist=False)
plt.savefig('Drug_Density_plot.svg')

# plt.figure()
# p2=sns.distplot(df[df['type'] == 'E']['number'], color=colors[0], label ='Drogas cuyos targets son genes causantes de endometriosis', hist=False)
# plt.savefig('Density_plot2.svg')
