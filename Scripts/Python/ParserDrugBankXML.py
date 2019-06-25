from bs4 import BeautifulSoup
from pprint import pprint
import xml.etree.cElementTree as etree
import sys

infile = open("full_database_30_04_2018.xml", "r")
# infile = open("aux/new.xml", "r")

parser = BeautifulSoup(infile.read(), "lxml")

print("Fichero leido")

drugs = parser.find_all("drug",created=True)

print("Fichero parseado");


output=open("drugbank.csv", "w")

output.write("#ID\tname\talias\tatcCodes\tsmiles\ttransporters\ttargets\tenzymes\tcarriers\tgroups\n");

atcCodes = {}

for drug in drugs:

	name=drug.find("name").string
	pprint(name)
	
	id=drug.find("drugbank-id", primary=True)
	smiles = ""
	transporters = []
	targets = []
	carriers = []
	enzymes = []
	syns = []
	atc = []
	groups = []


	cp = drug.select("calculated-properties > property")

	for p in cp:
		if p.kind.string == "SMILES":
			smiles = p.value.string
			break

	transporters_elems = drug.find_all("transporter")

	for t_elem in transporters_elems:
		aux = t_elem.find("gene-name")
		if aux and aux.string != None:
			transporters.append(aux.string)
	

	targets_elems = drug.find_all("target")

	for t_elem in targets_elems:
		aux = t_elem.find("gene-name")
		if aux and aux.string != None:
			targets.append(aux.string)


	groups_elems = drug.find_all("group")

	for group_elem in groups_elems:
		groups.append(group_elem.string)

	for t_elem in targets_elems:
		aux = t_elem.find("gene-name")
		if aux and aux.string != None:
			targets.append(aux.string)			

	enzymes_elems = drug.find_all("enzyme")


	for e_elem in enzymes_elems:
		aux = e_elem.find("gene-name")
		if aux and aux.string != None:
			enzymes.append(aux.string)

	carriers_elems = drug.find_all("carrier")

	for c_elem in carriers_elems:
		aux = c_elem.find("gene-name")
		if aux and aux.string != None:
			carriers.append(aux.string)

	synonyms_elems = drug.find_all("synonym")
	
	for syn_elem in synonyms_elems:
		syns.append(syn_elem.string)

	atc_codes_elems = drug.find_all("atc-code")

	for atc_elem in atc_codes_elems:
		code = atc_elem["code"]
		level_elems = atc_elem.find_all("level")
		atcCodes[code] = "-"

		atc_aux = [code]

		for level_elem in level_elems:
			code = level_elem["code"]
			value = level_elem.string
			atcCodes[code] = value
			atc_aux.append(code)

		atc.append(",".join(atc_aux))

	try:
		output.write(
			id.string.replace("\t", "").encode("utf-8").strip()  + "\t" + 
			name.replace("\t", "").encode("utf-8").strip() + "\t" + 
			";".join(set(syns)).replace("\t", "").encode('utf-8').strip()  + "\t" + 
			";".join(atc).replace("\t", "").encode('utf-8').strip()  + "\t" + 
			smiles.replace("\t", "").encode("utf-8").strip() + "\t"  + 
			";".join(set(transporters)).replace("\t", "").encode('utf-8').strip()  + "\t" + 
			";".join(set(targets)).replace("\t", "").encode('utf-8').strip()  + "\t" + 
			";".join(set(enzymes)).replace("\t", "").encode('utf-8').strip()  + "\t" + 
			";".join(set(carriers)).replace("\t", "").encode('utf-8').strip() + "\t" +
			";".join(set(groups)).encode('utf-8').strip() +
			"\n")

	except TypeError as e:
		# print drug
		pprint(id)
		pprint(smiles)
		pprint(transporters)
		pprint(targets)		
		pprint(groups)
		raise

	print("-----------")


atc_file = open("atc_codes.tsv", "w")

for key, value in atcCodes.iteritems():
	atc_file.write(key.replace("\t", "").encode("utf-8").strip() + "\t" + value.replace("\t", "").encode("utf-8").strip() + "\n")


atc_file.close()	
infile.close()
output.close()