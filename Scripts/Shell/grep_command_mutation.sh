 grep CM1417158 hgmd_pro_2019.1_hg38_annot_filter.vcf  | cut -f8 | sed 's/;/\n/g' | grep "^ANN" | cut -d"=" -f2 | sed 's/,/\n/g'

