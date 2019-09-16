#!/usr/bin/env bash
curl http://rest.kegg.jp/list/pathway/hsa  > pathways.csv

while read -r -a columns; do                   # read a line into an array
curl http://rest.kegg.jp/get/${columns[0]#*:}/kgml > ./kgml/${columns[0]#*:}

done < pathways.csv  
