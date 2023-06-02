For documentation only - files missing here:

"DHS_Index_and_Vocabulary_hg38_WM20190703.core.bed.gz"
"DHS_Index_and_Vocabulary_hg38_WM20190703.core.bed.gz.tbi"

Are not included due to file size, but can be created in the following way:
```
wget https://zenodo.org/record/3838751/files/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz
zcat DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz | tail -n+2 | sed 's/ /_/g' | awk '{OFS="\t"}{print $1,$8-1,$9,$4,$5,$10}' | sort -k1,2 -V | grep -v NA | bgzip > DHS_Index_and_Vocabulary_hg38_WM20190703.core.bed.gz 
tabix -p bed DHS_Index_and_Vocabulary_hg38_WM20190703.core.bed.gz 
```

".no_upload/Basic_QC_OK.txt" 
Has the following format:
s	is_female
Ind_ID1	TRUE
Ind_ID2	FALSE
Ind_ID3	FALSE
Ind_ID4	FALSE

".no_upload/covid_individuals.txt" 
Has the following format, where higher numbers in the column is_case corresponds to a more severe phenotype:
s	is_case	is_female
Ind_ID1	2	TRUE
Ind_ID2	1	FALSE
Ind_ID3	1	FALSE
Ind_ID4	1	FALSE

".no_upload/EURs_unrel.tsv" 
Has the following format, where higher numbers in the column is_case corresponds to a more severe phenotype, Age is self-explainatory, A1/B1/C1 are phenotype definitions:
s	is_case	is_female	Age	A1	B1	C1
Ind_ID1	2	TRUE	70	1	1	1
Ind_ID2	1	FALSE	70	0	1	1
Ind_ID3	1	FALSE	75	0	1	1
Ind_ID4	1	FALSE	75	0	1	1


