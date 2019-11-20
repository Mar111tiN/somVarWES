###Filter all exomes based on Kenichi criteria

#ExonicFunc.refGene != "unknown"
#misRate_tumor > 0.05
#P-value(EBCall) > 4
#P-value(fisher_realignment) > 1.3
#P-value(fisher) > 1.3
#misRate_normal < 0.02 (carefully confirm known driver mutations were not included)
#Also, germline SNPS with minor allele frequency of >=0.001 were removed, except for mutations with no calls in normal sample.

library(tidyverse)

#load raw Exome file###
Exome_Mut_all <- readr::read_tsv("~/Desktop/Japan/Bioinformatics/Analysis/PMBL Exomes/New Exomes/Raw_Files/201902_All_WES_except PMBL_44/Exome_mut_all.txt" , col_names = T, cols(Chr="c"))

#load sample filter files#
Exomes_passed_QC <-readr::read_tsv("~/Desktop/Japan/Bioinformatics/Analysis/Filter/Sample_Filter/Exomes_Passed_QC.txt", col_names = T)
Exomes_low_VAF <- readr::read_tsv("~/Desktop/Japan/Bioinformatics/Analysis/Filter/Sample_Filter/Exomes_below_VAF_threshold.txt" , col_names = T)


### Filter samples that passed QC and have median VAF >0.1####
WES_sample_filt <- Exome_Mut_all %>% 
  dplyr::filter(file_name %in% Exomes_passed_QC$file_name) %>% 
  dplyr::filter(!file_name %in% Exomes_low_VAF$file_name) %>% 
  dplyr::filter(!file_name == "I12_first")

#Filter otu unwanted regions#
WES_mut_filt_1 <- WES_sample_filt %>%
  dplyr::filter(ExonicFunc.refGene!="unknown" | is.na(ExonicFunc.refGene)) %>%
  dplyr::filter(!`AAChange.refGene`=="UNKNOWN" | is.na(`AAChange.refGene`)) %>%
  dplyr::filter(!Func.refGene %in% c("downstream","intergenic","intronic", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", 
                                     "ncRNA_splicing", "upstream", "upstream;downstream", "UTR3", "UTR5", "UTR5;UTR3")) 

#BASIC FILTERING CRITERIA#####
# A. moderate as a basic
WES_mut_filt_all_criteria <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3) %>% 
  dplyr::filter(depth_tumor >=10) %>% 
  dplyr::filter(`P-value(EBCall)`>= 4) %>% 
  dplyr::filter(`strandRatio_tumor` != 0 | `strandRatio_tumor` != 1) %>% 
  dplyr::filter(`misRate_tumor`>= 0.05) %>% 
  dplyr::filter(misRate_normal<= 0.05) %>% 
  dplyr::filter(`P-value(fisher)`>= 0.8 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 1.1 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

#If wanted Filter for SNPS that pop up in all databases
SNPs_annotated_in_any <- WES_mut_filt_all_criteria %>% 
  dplyr::filter(!is.na(`1000g2014oct_all`) |
                  !is.na(`1000g2010nov_all`) |
                  !is.na(`1000g2014oct_eur`) |
                  !is.na(esp6500siv2_all) |
                  (!is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`) & `HGVD_20131010:Frequency(NA/(NA+NR))` != "---") |
                  (!is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`) & `HGVD_20160412:Frequency(NA/(NA+NR))` != "---") |
                  (!is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`) & `ExAC:Frequency(AC_Adj/AN_Adj)` != "---") |
                  (!is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) & `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` != "---")) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(SNPs_annotated_in_any, "SNPs_annotated_in_any.txt")

### filter possibilities for SNPs

#A. Positive selection
filter(as.numeric(`1000g2014oct_all`) >0.001 & !is.na(`1000g2014oct_all`) & variantNum_normal!=0)
#OR#
dplyr::filter((as.numeric(`1000g2014oct_all`) > 0.001 & variantNum_normal>0)) 

#Negative selection
dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`))



#Basics all and excluding SNPs according to Kenichi
#A. strict
WES_mut_filt_all_criteria_no_SNP_strict <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=5) %>% 
  dplyr::filter(depth_tumor >=15) %>% 
  dplyr::filter(`P-value(EBCall)`>= 5.5) %>% 
  dplyr::filter(`strandRatio_tumor` != 0 | `strandRatio_tumor` != 1) %>% 
  dplyr::filter(`misRate_tumor`>= 0.07) %>% 
  dplyr::filter(misRate_normal<= 0.02) %>% 
  dplyr::filter(`P-value(fisher)`>= 2 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 2 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 


readr::write_tsv(WES_mut_filt_all_criteria_no_SNP_strict, "WES_filtered_no_SNP_strict.txt")


#B. moderate
WES_mut_filt_all_criteria_no_SNP_moderate <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3) %>% 
  dplyr::filter(depth_tumor >=10) %>% 
  dplyr::filter(`P-value(EBCall)`>= 4) %>% 
  dplyr::filter(`strandRatio_tumor` != 0 | `strandRatio_tumor` != 1) %>% 
  dplyr::filter(`misRate_tumor`>= 0.05) %>% 
  dplyr::filter(misRate_normal<= 0.05) %>% 
  dplyr::filter(`P-value(fisher)`>= 0.8 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 1.1 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 


readr::write_tsv(WES_mut_filt_all_criteria_no_SNP_moderate, "WES_filtered_no_SNP_moderate.txt")


#C. loose
WES_mut_filt_all_criteria_no_SNP_loose <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=2) %>% 
  dplyr::filter(depth_tumor >=8) %>% 
  dplyr::filter(`P-value(EBCall)`>=3) %>% 
  dplyr::filter(`misRate_tumor`>= 0.03) %>% 
  dplyr::filter(misRate_normal<= 0.07) %>% 
  dplyr::filter(`P-value(fisher)`>= 0.5 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 0.5 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) 


readr::write_tsv(WES_mut_filt_all_criteria_no_SNP_loose, "WES_filtered_no_SNP_loose.txt")

#excluded SNPs according to Kenichi
SNPs <- WES_mut_filt_all_criteria %>% 
  dplyr::filter((as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal>0) | 
                  (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal>0) |
                  (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal>0) |
                  (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal>0) |
                  (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal>0) |
                  (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal>0) |
                  (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal>0) |
                  (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal>0)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(SNPs, "SNPs.txt")

#pEB 3-4
pEB_3_to_4 <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3) %>% 
  dplyr::filter(depth_tumor >=10) %>% 
  dplyr::filter(`P-value(EBCall)`>= 3 & `P-value(EBCall)`<4) %>% 
  dplyr::filter(`strandRatio_tumor` != 0 | `strandRatio_tumor` != 1) %>% 
  dplyr::filter(`misRate_tumor`>= 0.05) %>% 
  dplyr::filter(misRate_normal<= 0.05) %>% 
  dplyr::filter(`P-value(fisher)`>= 0.8 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 1.1 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(pEB_3_to_4, "pEB_3_to_4.txt")

#P-FISHER <0.8
pFisher_smaller_0.8 <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3) %>% 
  dplyr::filter(depth_tumor >=10) %>% 
  dplyr::filter(`P-value(EBCall)`>= 4) %>% 
  dplyr::filter(`strandRatio_tumor` != 0 | `strandRatio_tumor` != 1) %>% 
  dplyr::filter(`misRate_tumor`>= 0.05) %>% 
  dplyr::filter(misRate_normal<= 0.05) %>% 
  dplyr::filter(`P-value(fisher)`>= 0 & `P-value(fisher)`<0.8 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 1.1 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(pFisher_smaller_0.8, "pFisher_smaller_0.8.txt")

#P-real 0.5-1.1####
P_real_0.5_to_1.1 <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3) %>% 
  dplyr::filter(depth_tumor >=10) %>% 
  dplyr::filter(`P-value(EBCall)`>= 4) %>% 
  dplyr::filter(`strandRatio_tumor` != 0 | `strandRatio_tumor` != 1) %>% 
  dplyr::filter(`misRate_tumor`>= 0.05) %>% 
  dplyr::filter(misRate_normal<= 0.05) %>% 
  dplyr::filter(`P-value(fisher)`>= 0.8 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >=0.5 & `P-value(fisher_realignment)` < 1.1 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(P_real_0.5_to_1.1, "P_real_0.5_to_1.1.txt")

unique(WES_mut_filt_1$file_name)

#VAF0.02-0.05
misRate_0.2_to_0.5 <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3) %>% 
  dplyr::filter(depth_tumor >=10) %>% 
  dplyr::filter(`P-value(EBCall)`>= 4) %>% 
  dplyr::filter(`strandRatio_tumor` != 0 | `strandRatio_tumor` != 1) %>% 
  dplyr::filter(`misRate_tumor`>= 0.02 & `misRate_tumor` < 0.05) %>% 
  dplyr::filter(misRate_normal<= 0.05) %>% 
  dplyr::filter(`P-value(fisher)`>= 0.8 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 1.1 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(misRate_0.2_to_0.5, "misRate_0.2_to_0.5.txt")

#misRate_normal 0.05-0.07
misRate_N_0.05_to_0.07 <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3) %>% 
  dplyr::filter(depth_tumor >=10) %>% 
  dplyr::filter(`P-value(EBCall)`>= 4) %>% 
  dplyr::filter(`strandRatio_tumor` != 0 | `strandRatio_tumor` != 1) %>% 
  dplyr::filter(`misRate_tumor` >= 0.05) %>% 
  dplyr::filter(misRate_normal <= 0.07 & misRate_normal > 0.05) %>% 
  dplyr::filter(`P-value(fisher)`>= 0.8 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 1.1 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(misRate_N_0.05_to_0.07, "misRate_N_0.05_to_0.07.txt")

#strand ratio 0 or 1
strand_ratio_0_1 <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3) %>% 
  dplyr::filter(depth_tumor >=10) %>% 
  dplyr::filter(`P-value(EBCall)`>= 4) %>% 
  dplyr::filter(`strandRatio_tumor` == 0 | `strandRatio_tumor` == 1) %>% 
  dplyr::filter(`misRate_tumor` >= 0.05) %>% 
  dplyr::filter(misRate_normal <= 0.05) %>% 
  dplyr::filter(`P-value(fisher)`>= 0.8 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 1.1 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`)) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(strand_ratio_0_1, "strand_ratio_0_1.txt")


#Rest of mutations
WES_Rest_mut <- WES_mut_filt_1 %>% 
  dplyr::filter(!Unique_Identifier %in% WES_mut_filt_all_criteria_no_SNP_moderate$Unique_Identifier) %>% 
  dplyr::filter(!Unique_Identifier %in% SNPs$Unique_Identifier) %>% 
  dplyr::filter(!Unique_Identifier %in% pEB_3_to_4$Unique_Identifier) %>% 
  dplyr::filter(!Unique_Identifier %in% misRate_0.2_to_0.5$Unique_Identifier) %>% 
  dplyr::filter(!Unique_Identifier %in% misRate_N_0.05_to_0.07$Unique_Identifier) %>% 
  dplyr::filter(!Unique_Identifier %in% pFisher_smaller_0.8$Unique_Identifier) %>% 
  dplyr::filter(!Unique_Identifier %in% P_real_0.5_to_1.1$Unique_Identifier) %>% 
  dplyr::filter (!Unique_Identifier %in% strand_ratio_0_1$Unique_Identifier) %>% 
  dplyr::add_count(Gene.refGene, name = "Frq_Gene", sort=TRUE) 

readr::write_tsv(WES_Rest_mut, "WES_Rest_mut.txt")







##### filter true somatic mutations for TP validation#####
WES_mut_filt_2 <- WES_mut_filt_1 %>% 
  mutate("VAF_T/N" = misRate_tumor/misRate_normal) 

True_som_for_TP <- WES_mut_filt_2 %>% 
  dplyr::filter(depth_tumor >10 & depth_normal >10) %>% 
  dplyr::filter(variantNum_tumor >=2) %>% 
  dplyr::filter(`misRate_tumor`>= 0.1) %>% 
  dplyr::filter(`VAF_T/N` >=5) %>% 
  dplyr::filter(Gene.refGene %in% Target_Genes$Gene.refGene)

Germline_SNP_for_TP <- WES_mut_filt_2 %>% 
  dplyr::filter(depth_tumor >10 & depth_normal >10) %>% 
  dplyr::filter(misRate_normal >= 0.15) %>% 
  dplyr::filter(Gene.refGene %in% Target_Genes$Gene.refGene)

Error1_for_TP <- WES_mut_filt_2 %>% 
  dplyr::filter(depth_tumor >10 & depth_normal >10) %>% 
  dplyr::filter(depth_tumor > 100) %>% 
  dplyr::filter(Gene.refGene %in% Target_Genes$Gene.refGene)

Error2_for_TP <- WES_mut_filt_2 %>% 
  dplyr::filter(depth_tumor >10 & depth_normal >10) %>% 
  dplyr::filter(variantNum_normal >=2) %>% 
  dplyr::filter(misRate_normal < 0.15) %>% 
  dplyr::filter(`VAF_T/N` < 5) %>% 
  dplyr::filter(Gene.refGene %in% Target_Genes$Gene.refGene)

Unambiuous_for_TP <- WES_mut_filt_2 %>% 
  dplyr::filter(Gene.refGene %in% Target_Genes$Gene.refGene) %>% 
  dplyr::filter(!Unique_identifier %in% True_som_for_TP$Unique_Identifier) %>% 
  dplyr::filter(!Unique_identifier %in% Germline_SNP_for_TP$Unique_Identifier) %>% 
  dplyr::filter(!Unique_identifier %in% Error1_for_TP$Unique_Identifier) %>% 
  dplyr::filter(!Unique_identifier %in% Error2_for_TP$Unique_Identifier) %>% 
  
  
  
  
  
  
  
  
  VNum3_Depth10 <- WES_mut_filt_1 %>% 
  dplyr::filter(variantNum_tumor >=3 & variantNum_tumor < 4 | (depth_tumor >=10 & depth_tumor <12)) %>% 
  dplyr::filter(`P-value(EBCall)`>= 4) %>% 
  dplyr::filter(`P-value(fisher)`>= 1.3 | is.na(`P-value(fisher)`)) %>% 
  dplyr::filter(`P-value(fisher_realignment)` >= 1.3 | is.na(`P-value(fisher_realignment)`)) %>% 
  dplyr::filter(`misRate_tumor`>= 0.05) %>% 
  dplyr::filter(misRate_normal<= 0.02) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_all`) <0.001 | (as.numeric(`1000g2014oct_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2010nov_all`) <0.001 | (as.numeric(`1000g2010nov_all`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2010nov_all`)) %>% 
  dplyr::filter(as.numeric(`1000g2014oct_eur`) <0.001 | (as.numeric(`1000g2014oct_eur`) >= 0.001 & variantNum_normal==0) | is.na(`1000g2014oct_eur`)) %>% 
  dplyr::filter(as.numeric(esp6500siv2_all) <0.001 | (as.numeric(esp6500siv2_all) >= 0.001 & variantNum_normal==0) | is.na(esp6500siv2_all)) %>% 
  dplyr::filter(as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20131010:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20131010:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20131010:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) <0.001 | (as.numeric(`HGVD_20160412:Frequency(NA/(NA+NR))`) >= 0.001 & variantNum_normal==0) | `HGVD_20160412:Frequency(NA/(NA+NR))` =="---" | is.na(`HGVD_20160412:Frequency(NA/(NA+NR))`)) %>% 
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_Adj/AN_Adj)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_Adj/AN_Adj)` =="---" | is.na(`ExAC:Frequency(AC_Adj/AN_Adj)`))  %>%
  dplyr::filter(as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) <0.001 | (as.numeric(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`) >= 0.001 & variantNum_normal==0) | `ExAC:Frequency(AC_POPMAX/AN_POPMAX)` =="---" | is.na(`ExAC:Frequency(AC_POPMAX/AN_POPMAX)`))

readr::write_tsv(VNum3_Depth10, "VNum3_Depth10.txt")
