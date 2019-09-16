# TFM
En la carpeta *Data*, se encuentra:  
**Información sobre GWAS y genes de enfermedad:** (carpeta DiseaseGenes)
- *Fichero con polimorfismos (SNPs) estrechamente asociados con la endometriosis que cumplen con el criterio (P < 1 × 10−5)*: GWAS_Associations.ods
- *Fichero con 150 genes de enfermedad (disease genes)*: GWAS_Genes.ods
- *Fichero con genes de enfermedad mapeados en KEGG*: Disease_Genes_mapped_in_KEGGPathways.ods

**Información sobre farmacos:** (carpeta Drugs)
- *Fichero con farmacos relacionados con endometriosis (seed drugs)*: Seed_drugs.ods
- *Fichero con farmacos que interactuan con los fármacos de endometriosis (seed + interacting drugs + disease gene-related drugs)*: DrugBank_drugs.ods
- *Fichero con la distribución de grado de los fármacos en la red final y en la red directa TD*: Degree_distribution_of_drugs.ods

**Información sobre genes diana (target):** (carpeta TargetGenes)
- *Fichero con 1502 genes diana (target genes)*: Drugbank_Target_Genes.ods
- *Fichero con genes diana mapeados en KEGG*: Target_Genes_mapped_in_KEGGPathways.ods

**Información sobre el interactoma de la base de datos KEGG:** (carpeta KEGG Interactome)
- *Fichero con los nodos y atributos del interactoma*: ATTS_genes.csv
- *Fichero con las relaciones y atributos del interactoma*: SIF_genes.csv
- *Fichero con todos los dataframes correspondientes a objetos KEGG (Entry, Relation, Pathway, Reaction)*: en carpeta dataframes

**Información sobre las variantes extraidas de la base de datos HGMD:** (carpeta HGMDVariants)
- *Fichero con las variantes y sus atributos de los genes diana de los fármacos priorizados*: Variants_filtered_by_Phenotype.ods

En la carpeta *Networks*, se encuentra una subcarpeta por cada subgrafo generado que contiene:

**Información sobre el subgrafo generado con *K* vecinos, con K=0..5:** (carpeta *K*Neighbors)
- *Fichero con los nodos y atributos del subgrafo*: ATTS_*K*Neighbors.csv
- *Fichero con las relaciones y atributos del subgrafo*: ATTS_*K*Neighbors.csv

**Información sobre los archivos de visualización de Cytoscape** (carpeta Cytoscape)
- *Red contexto molecular de la endometriosis*: All_Neighbors_Networks.cys
- *Red contexto farmacológico de la endometriosis*: Drug_network.cys
- *Red con las variantes de los genes diana de los fármacos priorizados*: Variant_network.cys
- *Red con los pathways del contexto molecular de la endometriosis*: Pathways.cys

En la carpeta *Scripts*, se encuentran todos los programas realizados para obtener los resultados del proyecto.  

Los scripts han sido realizados con Python 3.7.0 y R version 3.2.3.
