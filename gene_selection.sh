for i in `seq 1 5`
 do cut -f14 $biodb_human_cancer_gene_census_MapToEns75g_hg19_psl | sortu | while read chr
 do awk -v chr="$chr" '$14==chr' $biodb_human_cancer_gene_census_MapToEns75g_hg19_psl | shuf | h1
done
done  | sortu | cut -f10,14,16,17 | lest