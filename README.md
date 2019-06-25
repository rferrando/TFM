# TFM
En la carpeta *Data*, se encuentra:  
**Información sobre GWAS y genes de enfermedad:** (carpeta Disease Genes)
- *Fichero con 150 polimorfismos (SNPs) estrechamente asociados con la endometriosis que cumplen con el criterio (P < 1 × 10−5)*: GWAS_Info.csv
- *Fichero con 150 genes de enfermedad (disease genes)*: GWAS_by_gene_formatted.csv  

**Información sobre farmacos y genes diana (Target):** (carpeta Target Genes y Drugs)
- *Fichero con farmacos relacionados con endometriosis (seed drugs)*: seed_drugs_targets_endometriosis.csv
- *Fichero con farmacos que interactuan con los fármacos de endometriosis (seed + interacting drugs)*: approved_drugs_endometriosis_reduced.csv
- *Fichero con 1502 genes diana (target genes)*: Drugs_by_gene_TA.csv
- *Fichero con todas los fármacos aprobados en DrugBank*: drugbank_output.csv
- *Fichero con todas los fármacos aprobados que tienen como diana a algún gen de enfermedad.*: drugs_E.csv

**Información sobre el interactoma de la base de datos KEGG:** (carpeta KEGG Interactome)
- *Fichero con los nodos y atributos del interactoma*: ATTS_genes.csv
- *Fichero con las relaciones y atributos del interactoma*: ATTS_genes.csv

**Información sobre las variantes extraidas de la base de datos HGMD:** (carpeta HGMD Variants)
- *Fichero con las variantes y sus atributos para los genes de interés (Disease y Target)*: XXXXXX.csv

En la carpeta *Networks*, se encuentra una subcarpeta por cada subgrafo generado que contiene:
**Información sobre el subgrafo generado con *K* vecinos, con K=0..5:** (carpeta *K*Neighbors)
- *Fichero con los nodos y atributos del subgrafo*: ATTS_*K*Neighbors.csv
- *Fichero con las relaciones y atributos del subgrafo*: ATTS_*K*Neighbors.csv

En la carpeta *Scripts*, se encuentran todos los programas realizados para obtener los resultados del proyecto.  

Los scripts han sido realizados con Python 3.7.0 y R version 3.2.3.
