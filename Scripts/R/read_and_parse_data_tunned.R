# Script to parse KEGG database


# Settings ----------------------------------------------------------------

rm(list=ls())

# Set working directory
 setwd("/home/rosa/Master_Bioinformatica/TFM/KEGG/RData")

# install.packages("XML")
# Libraries
library(XML)

# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
library(biomaRt)

# Read and parse KGML data ------------------------------------------------

# Getting KGML files paths 
mypath = "/home/rosa/Master_Bioinformatica/TFM/KEGG/pathway"
myfiles = dir("/home/rosa/Master_Bioinformatica/TFM/KEGG/pathway")
length(myfiles)
# [1] 325

# Extract pathway KGML information

system.time(
  
  pathwaysInformation <- lapply(myfiles, function(xx){
    
    print(xx)
    
    x = paste(mypath,"/",xx,sep="")
    
    # Read file
    xmlpath = xmlParse(x)
    
    # First table: pathway information --
    
    # Get root
    root = xmlRoot(xmlpath)
    
    # Get information
    name = gsub("path:","", xmlGetAttr(xmlRoot(xmlpath),"name"),fixed=T)
    title = xmlGetAttr(xmlRoot(xmlpath),"title")
    link = xmlGetAttr(xmlRoot(xmlpath),"link")
    image = xmlGetAttr(xmlRoot(xmlpath),"image")
    
    # # Add information about group and subgroup
    # group = pathInfo$Group[pathInfo$pathwayID==gsub(".xml","",xx,fixed=T)]
    # subgroup = pathInfo$Subgroup[pathInfo$pathwayID==gsub(".xml","",xx,fixed=T)]
    
    # pathway table
    # pathway = data.frame(name=name,title=title,link=link,png=image,group=group,subgroup=subgroup)
    pathway = data.frame(name=name,title=title,link=link,png=image, stringsAsFactors = F)
    
    # Second table: entry information --
    
    en = getNodeSet(xmlpath,"//pathway/entry")
    
    entries = do.call("rbind",lapply(1:length(en),function(y){
      #     print(y)
      id = xmlGetAttr(en[[y]], "id")
      name = xmlGetAttr(en[[y]], "name") 
      type = xmlGetAttr(en[[y]], "type")
      reaction = xmlGetAttr(en[[y]],"reaction")
      if(is.null(reaction)){
        reaction="--"
      }
      link = xmlGetAttr(en[[y]], "link")
      if (is.null(link)){
        link = "--"
      }
      grap = getNodeSet(en[[y]], "graphics")   
      graphics_name = xmlGetAttr(grap[[1]],"name")
      if (is.null(graphics_name)){
        graphics_name = "--"
      }
      xx = xmlGetAttr(grap[[1]],"x")
      if (is.null(xx)){
        xx = "--"
      }
      yy = xmlGetAttr(grap[[1]], "y")
      if (is.null(yy)){
        yy="--"
      }
      if (type == "group"){
        comps = getNodeSet(en[[y]],"component")
        name = paste(unlist(lapply(1:length(comps),function(z){
          xmlGetAttr(comps[[z]],"id")
        })),collapse="&")
      }
      
      data.frame(id=id, genes=name ,type=type, reaction=reaction, link=link, 
                 name=graphics_name, x=xx, y=yy, stringsAsFactors = F)
      
    }))
    
    # Third table: relations information --
    
    rels = getNodeSet(xmlpath,"//pathway/relation")
    if(length(rels)==0){
      relations = NULL
    }else{
      relations = do.call("rbind",lapply(1:length(rels),function(y){
        #       print(y)
        entry1 = xmlGetAttr(rels[[y]],"entry1")
        entry2 = xmlGetAttr(rels[[y]], "entry2")
        type = xmlGetAttr(rels[[y]], "type")
        
        subt = getNodeSet(rels[[y]],"subtype")
        if (length(subt) == 0){
          subtype_name = subtype_value = "--"
        }else{
          subtype_name = paste(unlist(lapply(1:length(subt), function(z){
            xmlGetAttr(subt[[z]],"name")
          })),collapse=" ")
          subtype_value = paste(unlist(lapply(1:length(subt), function(z){
            xmlGetAttr(subt[[z]],"value")
          })),collapse=" ")
        }
        
        data.frame(entry1=entry1,entry2=entry2,type=type,subtype=subtype_name,
                   value=subtype_value, stringsAsFactors = F)
      }))
    }
    
    # # If there are groups
    # if (!is.null(relations) & any(entries$type=="group")){
    #   for (i in which(entries$type=="group")){
    #     # Add relations related to groups
    #     genes = unlist(strsplit(as.character(entries$genes[i]), split=" "))
    #     combs = combn(genes,2)
    #     new_relations_1 = do.call("rbind", lapply(1:ncol(combs), function(z){
    #       data.frame(entry1 = combs[1,z],
    #                  entry2 = combs[2,z],
    #                  type = "added",
    #                  subtype = "binding/association",
    #                  value = "-new-",
    #                  stringsAsFactors = F)
    #     }))
    #     new_relations_2 = do.call("rbind", lapply(genes, function(z){
    #       data.frame(entry1 = z,
    #                  entry2 = entries$id[i],
    #                  type = "added",
    #                  subtype = "scaffold",
    #                  value = "-new->",
    #                  stringsAsFactors = F)
    #     }))
    #     relations = rbind(relations, new_relations_1, new_relations_2)
    #     # Change genes column in entries
    #     entries$genes[i] = paste(entries$genes[entries$id%in%genes], collapse = " ")
    #     # Change name
    #     new_name = unlist(lapply(genes, function(z){
    #       unlist(strsplit(entries$name[entries$id==z], split = ", "))[1]
    #     }))
    #     entries$name[i] = paste(new_name, collapse = " ")
    #   }
    # }
    # 
    
    # Fourth table ###
    
    reac = getNodeSet(xmlpath,"//pathway/reaction")
    if (length(reac) == 0){
      reactions = NULL
    }else{
      reactions = do.call("rbind",lapply(1:length(reac),function(y){
        id = xmlGetAttr(reac[[y]],"id")
        name = xmlGetAttr(reac[[y]], "name")
        type = xmlGetAttr(reac[[y]], "type")
        subs = getNodeSet(reac[[y]],"substrate")
        subs_id = paste(unlist(lapply(1:length(subs),function(z){
          xmlGetAttr(subs[[z]],"id")
        })),collapse=" + ")
        prod = getNodeSet(reac[[y]],"product")
        prod_id = paste(unlist(lapply(1:length(prod),function(z){
          xmlGetAttr(prod[[z]],"id")
        })),collapse=" + ")
        substrate_name= paste(unlist(lapply(1:length(subs),function(z){
          xmlGetAttr(subs[[z]],"name")
        })),collapse=" + ")
        product_name = paste(unlist(lapply(1:length(prod),function(z){
          xmlGetAttr(prod[[z]],"name")
        })),collapse=" + ")
        
        data.frame(enzyme = id, id = name, type = type, substrate=subs_id, subs2 = substrate_name, 
                   products = prod_id, prods2 = product_name, stringsAsFactors = F)
      }))
    }
    
    
    
    list(pathway = pathway,
         entries=entries,
         relations=relations,
         reactions=reactions)
  })
)

# user  system elapsed 
# 47.608   0.044  47.532 

# Names KEGG info
names(pathwaysInformation) = unlist(lapply(myfiles,function(x){
  gsub(".xml","",x)
}))

# Create table gene-pathway

system.time(
  
  annot_table <- do.call("rbind", lapply(pathwaysInformation, function(x){
    genes = unique(unlist(strsplit(as.character(x$entries$genes[x$entries$type=="gene"]),split=" ")))
    data.frame(pathway = x$pathway$name, 
               gene = genes, 
               stringsAsFactors = F)
  }))

)

# user  system elapsed 
# 0.160   0.004   0.165 

rownames(annot_table) = NULL

# Move to GeneName identifiers --------------------------------------------

# Set mart
ensembl=useMart("ensembl")

# Set dataset
human = useDataset("hsapiens_gene_ensembl",mart=ensembl)

entrez_genes = gsub("hsa:", "", unique(annot_table$gene), fixed=T)

system.time(
  
  genes_corres <- getBM(attributes=c('entrezgene', 'hgnc_symbol'), 
                       filters = 'entrezgene', 
                       values = entrez_genes, 
                       mart = human)

)
# user  system elapsed 
# 0.028   0.004   3.365 

# Remove genes with no annotation in hgnc_symbol
genes_corres = genes_corres[genes_corres$hgnc_symbol!="",]

# How many entrezgene are not annotated to geneName?
length(setdiff(entrez_genes, unique(genes_corres$entrezgene))) #40
length(setdiff(entrez_genes, genes_corres$entrezgene))*100/length(entrez_genes) # 0.5344735

# Annotate them manually - From NCBI gene database (entrezGene)
manual_annot = data.frame(rbind(c("1350", "COX7C"),
                                c("25885", "POLR1A"),
                                c("390877", "LOC390877"),
                                c("606495", "CYB5RL"),
                                c("113189", "CHST14"),
                                c("10775", "POP4"),
                                c("100288562", "LOC100288562"),
                                c("101929601", "LOC101929601"),
                                c("101929627", "LOC101929627"),
                                c("101930111", "LOC101930111"),
                                c("643802", "LOC643802"),
                                c("10607", "TBL3"),
                                c("6207", "RPS13"),
                                c("6202", "RPS8"),
                                c("100529097", "RPL36A-HNRNPH2"),
                                c("7515", "XRCC1"),
                                c("102723407", "LOC102723407"),
                                c("84876", "ORAI1"),
                                c("4985", "OPRD1"),
                                c("100132074", "FOXO6"),
                                c("2559", "GABRA6"),
                                c("100996746", "SPDYE11"),
                                c("102723849", "SPDYE17"),
                                c("3211", "HOXB1"),
                                c("567", "B2M"),
                                c("107181291", "POP3"),
                                c("4651", "MYO10"),
                                c("26492", "OR8G2P"),
                                c("26687", "OR4E1"),
                                c("246721", "POLR2J2"),
                                c("81691", "REXO5"),
                                c("102800317", "LOC400927-CSNK1E"),
                                c("100996758", "NPY4R2"),
                                c("105379861", "LOC105379861"),
                                c("102724652", "LOC102724652"),
                                c("102723996", "LOC102723996"),
                                c("102723532", "LOC102723532"),
                                c("102724428", "LOC102724428")), 
                          stringsAsFactors = F)
colnames(manual_annot) = colnames(genes_corres)

genes_corres = rbind(genes_corres, manual_annot)

# How many entrezGenes are annotated twice?
sum(duplicated(genes_corres$entrezgene)) # [1] 44
# --> There are 18 entrezGene ids that map to more than one geneName --> We decide to duplicate the information

annot_table_GN = do.call("rbind", lapply(unique(annot_table$pathway), function(x){
  genes = gsub("hsa:", "", annot_table$gene[annot_table$pathway==x], fixed=T)
  data.frame(pathway = x,
             gene = unique(genes_corres$hgnc_symbol[genes_corres$entrezgene%in%genes]),
             stringsAsFactors = F)
}))

# Table describing the variables
description = data.frame(objectName = c("pathwaysInformation", "annot_table", "genes_corres", "annot_table_GN"),
                         description = c("Information extracted from KGML files",
                                         "Annotation gene-pathway (geneID)",
                                         "Correspondence entrezGene - geneName",
                                         "Annotation gene-pathway (geneName)"),
                         class = c("list", "data.frame", "data.frame", "data.frame"),
                         stringsAsFactors = F)


save(pathwaysInformation, annot_table, genes_corres, annot_table_GN,
     file="KEGG_database_info_tunned.RData")
