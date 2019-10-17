# TFM
En la carpeta ***Data***, se encuentra:  
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
- *Fichero con los nodos y atributos de la info de pathways*: ATTS_pathways.csv
- *Fichero con las relaciones y atributos de la info de pathways*: SIF_pathways.csv
- *Fichero con todos los dataframes correspondientes a objetos KEGG (Entry, Relation, Pathway, Reaction)*: en carpeta dataframes

**Información sobre las variantes extraidas de la base de datos HGMD:** (carpeta HGMDVariants)
- *Fichero con las variantes y sus atributos de los genes diana de los fármacos priorizados*: Variants_filtered_by_Phenotype.ods
- *Fichero con los fenotipos de variantes priorizados*: Phen_priorizated.csv


En la carpeta ***Networks***, se encuentra una subcarpeta por cada subgrafo generado que contiene:  
**Información sobre el subgrafo generado con *K* vecinos, con K=0..5:** (carpeta *K*Neighbors)
- *Fichero con los nodos y atributos del subgrafo*: ATTS_*K*Neighbors.csv
- *Fichero con las relaciones y atributos del subgrafo*: ATTS_*K*Neighbors.csv

**Información sobre los archivos de visualización de Cytoscape** (carpeta Cytoscape)
- *Red contexto molecular de la endometriosis*: Red_final_2_vecinos.cys
- *Red contexto farmacológico de la endometriosis*: Drug_network.cys
- *Red con las variantes de los genes diana de los fármacos priorizados*: Variant_network.cys
- *Red con los pathways del contexto molecular de la endometriosis*: Pathways.cys

**Información sobre los filtros utilizados para la explotación del modelo** (carpeta FilterAnalysis)
- *Filtro para estudio de proximidad*: Disease_conectado_a_al_menos_2_TP y TS_conectados_a_al_menos_2_TP_y_al_menos_1_Disease
- *Filtro para estudio farmacológico*: Drogas_con_degree >= 35_en_final_y_>= 8_en_directa-TD

En la carpeta ***Scripts***, se encuentran:  
**Todos los programas realizados para obtener los resultados del proyecto**

Los scripts han sido realizados con programación Shell, Python 3.7.0 y R version 3.2.3.

- *Anotación de todas las variantes de HGMD*: variant_annot_snpEff.sh
- *Búsqueda de las variantes ligadas a nuestros genes de interés en HGMD*: AccessHGMDVariants.py
- *Filtrado de las variantes ligadas a nuestros genes de interés en el fichero VCF anotado*: variant_filter_bcftools.sh
- *Procesado(parsing) de la bbdd de DrugBank - Output: fichero .csv con formato ID DrugBank | Nombre del fármaco | targets | enzymes | transporters | carriers, donde cada fila representa un fármaco y las proteínas que interactúan con él*: ParserDrugBankXML.py
- *Extracción de todos los fármacos aprobados (seed drugs e interacting drugs) que mapearemos en nuestra red, indicando para cada uno de ellos el número (degree) de nodos target con los que interactúa*: ParseDrugBankInfo.py
- *Extracción de todos los genes target, transporter, enzyme y carrier de los fármacos aprobados, configurando la lista final de genes diana (target)*: ParseRelatedTargetsInfo.py
- *Extracción de todos los genes reportados y mapped
gene(s) de los catalogos GWAS, configurando la lista de genes disease*: ParseGWASInfo.py
- *Anotación de los genes disease y target, accediendo a los servicios web de MyGene.info para consultar con cada identificador gene symbol su identificador Entrez y Ensembl, teniendo en cuenta todos sus alias*: MapGeneSymbolToEntrez.py
- *Acceso a la base de datos KEGG PATHWAY a través de su API pública: descarga de los ficheros KGML*: download_kgml.sh
- *Procesado (parsing) de los ficheros KGML, obteniendo un dataframe por cada objeto: Pathway, Entry, Relation y Reaction*: read_and_parse_data_tunned.r
- *Obtención de los ficheros de nodos (ATTS) y relaciones (SIF) de la red de genes y de la red de pathways*: tranform_KEGG_data_to_network.r
- *Generación de los ficheros ATTS y SIF para cada uno de los subgrafos (n = 0,1,2,3,4,5)*: Create_Directed_SubNetwork_for_Clustering.py
- *Anotación los genes que estaban en nomenclatura KEGG (hsa+EntrezID) a gene symbol para usar esa anotación como etiquetas de los nodos en Cytoscape*: Convert_nodes_to_symbol.py
- *Cálculo y representación de los máximos relativos de las relaciones de nodos disease con primary targets y nodos secondary target con primary target. En el caso de la distribución de grado de los nodos disease con los primary target, el 94 % de los nodos disease no están conectados con ningún nodo primary target por lo que el máximo se ha calculado teniendo en cuenta sólo los que están directamente conectados*: TS_Analysis_Boxplot.py
- *Integración de los fármacos en la red final de N vecinos, añadiendo nuevos nodos a la red (seed drugs, interacting drugs y disease gene-related drugs) y añadiendo nuevas relaciones, representadas por las interacciones de esos fármacos con sus genes diana, (target y disease)*: drug_target_interaction_network.py y
adding_disease_drug_interaction_to_targets.py
- *Cálculo y representación de los máximos relativos de las relaciones de llos fármacos con sus genes diana*: Drug_Analysis_Boxplot.py
- *Obtención de la intersección de los fármacos priorizados en el estudio farmacológico con los fármacos que interactúan con las dianas terapeúticas potenciales priorizadas en el estudio de proximidad*: get_drugs_intersected.py
- *Obtención de las variantes de las dianas de los fármacos intersectados filtrando por los fenotipos priorizados*: Priorization_of_variants.py
- *Integración de las variantes en le red final de N vecinos con sus fármacos*: Adding_variants_to_targets.py
