
module load snpEff

snpeff="java -jar /soft/apps/snpEff/4.3s-Java-1.8.0_92/snpEff/snpEff.jar "

input=hgmd_pro_2019.1_hg38.vcf
output=hgmd_pro_2019.1_hg38_annot.vcf

$snpeff ann GRCh38.86 $input >  $output

module unload snpEff
