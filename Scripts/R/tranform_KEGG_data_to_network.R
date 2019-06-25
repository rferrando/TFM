########## USO DE LIBRERIAS ####################################

library(dplyr)
library(stringr)
library(plyr)
# install.packages("tidyverse")
library(tidyverse)

########## DIRECTORIO DE TRABAJO ####################################

# Set working directory
setwd("/home/rosa/Master_Bioinformatica/TFM/KEGG/")
dir=getwd()


########## IMPORT DATOS DE KEGG Y CREACION DE FICHEROS RDATA ####################################

# Cargo los datos de KEGG
load(paste(dir,"/RData/KEGG_database_info_tunned.RData",sep=""))

#Cargo los datos de los pathways clasificados por enfermedad o no enfermedad
load(paste(dir,"/RData/KEGG_annotation_January2019.RData",sep=""))

#Filtro los pathway que no son de enfermedad:
pathways_noDisease = pathwaysInfo[grep("FALSE", pathwaysInfo$disease),]

# modificar el pathways information original, quedandome solo con los pathways de la lista que no son enfermedad ------------------------------
#De este modo el resto de dataframes no incluirán los pathways de enfermedades
pathwaysInformation_noDisease = pathwaysInformation[pathways_noDisease$name]

# Construyo el dataframe de los pathways
pathways = do.call("rbind", lapply(pathwaysInformation_noDisease, FUN=function(pathway){pathway$pathway[c(1,2)]}))
save(pathways,file=paste(dir,"/RData/pathways.RData",sep=""))

# Construyo el dataframe de los entries
entries = do.call("rbind", lapply(pathwaysInformation_noDisease, FUN=function(pathway){pathway$entries[c(1,2,3)]}))
entries$pathway = lapply(rownames(entries), FUN=function(pathway){strsplit(pathway, ".", fixed=TRUE)[[1]][1]})
save(entries,file=paste(dir,"/RData/entries.RData",sep=""))

#Lo uso para poder filtrar por la columna pathway desde la consola R Studio
entries_test = entries
entries_test$pathway = unlist(entries_test$pathway, recursive = TRUE, use.names = TRUE)

# Construyo el dataframe de las relaciones
relations = do.call("rbind", lapply(pathwaysInformation_noDisease, FUN=function(pathway){pathway$relations}))
relations$pathway = lapply(rownames(relations), FUN=function(pathway){strsplit(pathway, ".", fixed=TRUE)[[1]][1]})

relations_test = relations
relations_test$pathway = unlist(relations_test$pathway, recursive = TRUE, use.names = TRUE)
save(relations,file=paste(dir,"/RData/relations.RData",sep=""))

# Construyo el dataframe de las reacciones
reactions = do.call("rbind", lapply(pathwaysInformation_noDisease, FUN=function(pathway){pathway$reactions}))
reactions$pathway = lapply(rownames(reactions), FUN=function(pathway){strsplit(pathway, ".", fixed=TRUE)[[1]][1]})
save(reactions,file=paste(dir,"/RData/reactions.RData",sep=""))

##########    CONSTRUCCIÓN RED DE GENES ####################################

# Construyo la red de genes, quitando las entries de tipo "map"
entries_genes_filtered <- entries %>% filter(type != "map")

#Eliminamos también las entries de tipo ortologo, compound y brite, para quedarnos sólo con las de tipo gen y grupo
entries_genes_filtered<- entries_genes_filtered %>% filter(type != "ortholog")
entries_genes_filtered<- entries_genes_filtered %>% filter(type != "compound")
entries_genes_filtered<- entries_genes_filtered %>% filter(type != "brite")

# Segrego en dos dataframe distintos las entries de tipo gen y las de tipo grupo
genes_group=entries_genes_filtered[entries_genes_filtered$type=="group",1:4] 
genes_gene=entries_genes_filtered[entries_genes_filtered$type=="gene",1:4] 

#A las entries de tipo grupo, les aplico una función para devolver la lista de genes de cada grupo en vez del ID (cada grupo separado por & - AND)
genes_group$list_genes = apply(genes_group, 1, FUN=function(group){paste(lapply(strsplit(group$genes, "&", fixed=TRUE)[[1]],FUN=function(genID){entries_genes_filtered[entries_genes_filtered$id == genID & entries_genes_filtered$pathway == group$pathway, c("genes")]}), collapse = "&")})
#Para las entries de tipo gen la lista de genes será la que venía en la columna genes (cada gen separado por espacio - OR)
genes_gene <- mutate(genes_gene, list_genes = genes)

#En el dataframe entries_genes combinaré todas las entries con sus listas de genes
entries_genes = rbind(genes_group, genes_gene)
entries_genes <- entries_genes[order(entries_genes$list_genes),]
entries_genes$pathway = unlist(entries_genes$pathway, recursive = TRUE, use.names = TRUE)

save(entries_genes,file=paste(dir,"/RData/entries_genes.RData",sep=""))
write.csv(entries_genes, file="entries_genes.csv",row.names=TRUE)

##########    FICHERO ATTS --> RED DE GENES ####################################

#Agrupa (group_by) y sumariza (summarise) cada entry (list_genes) con los pathways en los que interviene
counts_entries <- ddply(entries_genes, .(list_genes), summarize, count = paste(pathway, collapse = ";"))

counts_entries$pathways = apply(counts_entries, 1, FUN=function(entry){
   x = strsplit(as.character(entry["count"]), ";", fixed=TRUE)[[1]]
   paste(x[!duplicated(x)], collapse = ";")
  })

ATTS_genes = counts_entries[,c("list_genes", "pathways")]

save(ATTS_genes,file=paste(dir,"/RData/ATTS_genes.RData",sep=""))
write.csv(ATTS_genes, file="ATTS_genes.csv",row.names=TRUE)


##########    FICHERO SIF --> RED DE GENES ####################################

relations_genes = relations[,1:6]
dim(relations_genes)

#Mapearemos los ID de relations_genes con entries_genes para obtener la lista de genes por cada ID
relations_genes$list_genes1 = apply(relations_genes, 1, FUN=function(entry_gen){
  entries_genes[entries_genes$id==entry_gen$entry1 & entries_genes$pathway==entry_gen$pathway, c("list_genes")]
})
relations_genes$list_genes2 = apply(relations_genes, 1, FUN=function(entry_gen){
  entries_genes[entries_genes$id==entry_gen$entry2 & entries_genes$pathway==entry_gen$pathway, c("list_genes")]
})

#Me quedo solo con aquellas entries que estén en mi fichero ATTS_genes
relations_genes<- relations_genes %>% filter(lengths(list_genes1) > 0)
relations_genes<- relations_genes %>% filter(lengths(list_genes2) > 0)
dim(relations_genes)

#Para aquellas relaciones de subtipo Compound, voy a evaluar si se  trata de una ruta dirigida o bidireccional
#observando el tipo de las reacciones que hay a ambos lados del compound
#Relation: entry1 entry2 compound
#Reaction1: enzyme = entry1 / product = compound --> tipo: reversible/irreversible
#Reaction2: enzyme = entry2 / substrate = compound --> tipo: reversible/irreversible
#El tipo de la relación corresponderá al tipo de la reacción 2 en cualquier casuística (REV/IRREV, IRREV/IRREV, REV/REV, IRREV/REV)
for (row in 1:nrow(relations_genes)) {
  subtype <- relations_genes[row, "subtype"]
  # print(subtype)
  # if(subtype != "compound") {
  if (is.na(str_extract(subtype, regex("compound")))) {
    next
  } else {
    subtype = "compound"
    if (!is.na(str_split(subtype, "compound ")[[1]][2])) {
      relations_genes[row, c("subtype")] <- as.character(str_split(subtype, "compound ")[[1]][2])
      next
    }
  }
  pathway <- relations_genes[row, "pathway"]
  
  # print(pathway[[1]])
  entry1 <- relations_genes[row, "entry1"]
  # print(entry1)
  entry2 <- relations_genes[row, "entry2"]
  # print(entry2)
  compound <- relations_genes[row, "value"]
  # print(compound)
  
  df<- reactions %>% filter(pathway == pathway[[1]] , enzyme == entry1, compound %in% (str_extract_all(products, "\\d+")[[1]]))
  
  if (nrow(df) != 0){   relations_genes[row, c("Type_entry1")] <-df$type
  print(relations_genes[row, c("Type_entry1")])
  } else {relations_genes[row, c("Type_entry1")] <- "irreversible"}
  
  df2<- reactions %>% filter(pathway == pathway[[1]] , enzyme == entry2, compound %in% (str_extract_all(substrate, "\\d+")[[1]]))
  print(df2)
  if (nrow(df2) != 0){
    relations_genes[row, c("Type_entry2")] <-df2$type
    print(relations_genes[row, c("Type_entry2")])
  } else {relations_genes[row, c("Type_entry2")] <- "irreversible"}

  relations_genes[row, c("Type_compound")] <- as.character(relations_genes[row, c("Type_entry2")])
  relations_genes[row, "subtype"] <- paste(subtype, relations_genes[row, c("Type_compound")])
}

relations_genes <- apply(relations_genes,2,as.character)
write.csv(relations_genes, file="relations_genes.csv",row.names=TRUE)


SIF_genes = relations_genes[, c("list_genes1", "list_genes2", "subtype", "pathway")]
SIF_genes <- apply(SIF_genes,2,as.character)
write.csv(SIF_genes, file="SIF_genes.csv",row.names=TRUE)



##########    CONSTRUCCIÓN RED DE PATHWAYS ####################################

# Construyo la red de pathways, quedandome solo con las entries de tipo "map"
entries_pathways_filtered <- entries %>% filter(type == "map")
##Otra forma de eliminar en las entries los pathways relacionados con enfermedades (hsa:05 y 6.7 Endocrine and metabolic diseases) 
#a este nivel sería la siguientes (aunque como ya lo hacemos en pathwaysInformation de forma general ya no hace falta)
#diseases = c("path:hsa04930", "path:hsa04940", "path:hsa04950", "path:hsa04932", "path:hsa04931", "path:hsa04933", "path:hsa04934", "path:hsa01521", "path:hsa01523", "path:hsa01524")
#entries_path <- entries %>% filter(type == "map", !(str_detect(genes, "^path:hsa05")), !(genes %in% diseases))

entries_pathways <- entries_pathways_filtered[order(entries_pathways_filtered$genes),]
entries_pathways$pathway = unlist(entries_pathways$pathway, recursive = TRUE, use.names = TRUE)

save(entries_pathways,file=paste(dir,"/RData/entries_pathways.RData",sep=""))
write.csv(entries_pathways, file="entries_pathways.csv",row.names=TRUE)

##########    FICHERO ATTS --> RED DE PATHWAYS ####################################
ATTS_pathways = pathways

save(ATTS_pathways,file=paste(dir,"/RData/ATTS_pathways.RData",sep=""))
write.csv(ATTS_pathways, file="ATTS_pathways.csv",row.names=TRUE)


##########    FICHERO SIF --> RED DE PATHWAYS ####################################

SIF_pathways_duplicados = data.frame(pathway1 = entries_pathways$pathway,
                 pathway2 = entries_pathways$genes,
                 subtype = entries_pathways$type,
                 stringsAsFactors = F)

#En  SIF_pathways$pathway2 (que seran nuestros hsa:pathways) eliminamos el prefijo path para tener el pathway target
SIF_pathways_duplicados$pathway2 = gsub("path:", "", SIF_pathways_duplicados$pathway2)

# Eliminar duplicados en SIF ----------------------------------------------

#Eliminamos los ciclos
SIF_pathways_duplicados <- SIF_pathways_duplicados[SIF_pathways_duplicados$pathway1 != SIF_pathways_duplicados$pathway2, ]

#Eliminamos los duplicados tipo: glucolisis esta relacionado con b, c", pero en b y c tambien esta la glucolisis 
#Para ello creamos una nueva columna que tenga la union depathway1 con pathway2 de forma ordenada
#ejemplo: 
#1. Si aparece hsa20 con hsa10, te lo ordena con sort() >>> nueva columna= hsa10_hsa20
#1. Si aparece hsa10 con hsa20, te lo ordena con sort() >>> nueva columna= hsa10_hsa20
#de este modo ya puedes hacer un !duplicated segun esa nueva columna

SIF_pathways_duplicados$aux <- apply(SIF_pathways_duplicados, 1, function(row){ #creamos la nueva columna aplicando para cada fila la cond
  a <- row[["pathway1"]]
  b <- row[["pathway2"]]
  c <- sort(c(a,b))
  return (paste(c, collapse  = "_"))
})

SIF_pathways <- SIF_pathways_duplicados[!duplicated(SIF_pathways_duplicados$aux), ]  

SIF_pathways = SIF_pathways [,-4] 
write.csv(SIF_pathways, file="SIF_pathways.csv",row.names=TRUE)


# ####guardamos el global enviroment: -------------------------------------
save.image(file=paste(dir,"/RData/all_network_info.RData",sep=""))



