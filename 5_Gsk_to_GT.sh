### LT Dunning & JK Olofsson & D Wood  07/03/2017 ###
## Gskm_to_gene.sh v.1.0 ##
## This script maps data to a reference, and generates consensus sequences for the data (i.e. calls SNPs & non-variant positions)
## This script can be used for genome assemblies or for other types of reference (e.g. transcriptome), as long as the regions of interest are known 
## once you have added the path for each file then you should be able to run the all the way through 

## INFILES ##
## 1. bam files with data mapped to reference with SAMPLEID.bam (no underscores in sample ID; mapping quality filtered (Q>20); if using whole genome data it can speed up the process if you cut the bam file down to your regions of interest using bedtools intersect and a bed file of regions of interest (infile 2. below)).
BAMS=/home/$USER/map/1_proper_pairs

## 2. Bed file of regions of interest e.g. of BED file below for full gene windows in the genome:
#	Sc6ifeU_2050	184963	186318
#	Sc6ifeU_498	76624	81521
#	Sc6ifeU_498	100271	106291
#	Sc6ifeU_2420	632305	643083
# 	...
BED=merged_Oe6_covfilt_allcov_80per_mac3.bed

## 3. file listing all samples to be processed (N.B. one IDper line, no underscores in names, and above bam files need to be SAMPLEID.bam) 
SAMPLES=samples.txt

## 4. the reference fasta file used for mapping
REF=/home/$USER/olive_refgenome_Oe6

## 5. Python scrip  countBases2.py
countBases2=countBases2.py

## 6. Ref gene file (basically a bed file with an extra column at the start indicating which gene things belong to, if its a transcript you can just duplicate the first column)
# e.g. 
#	GeneA	Sc6ifeU_2050	184963	186318
#	GeneB	Sc6ifeU_498	100271	106291
#	GeneC	Sc6ifeU_2420	632305	643083
# 	..

############

## Step 1: make a gene file for each position you want to call from the BED file, resultant file should look like this (can take a little while... Optimised by Danny):
#	Sc6ifeU_2050	184963		
#	Sc6ifeU_2050	184964
#	Sc6ifeU_2050	184965
#	Sc6ifeU_2050	184966
# 	...
perl split_annotated_seq.pl ${BED}
## Step 2: make this file more useful for downstream work, e.g. 
#	Sc6ifeU_2050	184963	184964	/	+	
#	Sc6ifeU_2050	184964	184965	/	+
#	Sc6ifeU_2050	184965	184966	/	+
#	Sc6ifeU_2050	184966	184966	/	+
# 	...
mkdir manual_GT
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' Step_1.txt > manual_GT/Step_2.txt
rm Step_1.txt
awk 'NR%2==0' manual_GT/Step_2.txt > manual_GT/temp && mv manual_GT/temp manual_GT/Step_2.txt

# Step 3: make mpileup file using samtools and play around with the formatting so that every position has a call
cat ${SAMPLES} | while read line ; do samtools mpileup -A -l manual_GT/Step_2.txt -f ${REF}/Oe6.scaffolds.fa ${BAMS}/"$line"_Oe6_bowtie_unique_f2_sort.bam > manual_GT/"$line"_Oe6.mpileup ; done
cat ${SAMPLES} | while read line ; do awk 'BEGIN{OFS="\t"}; {print $1":"$2, $0}' manual_GT/"$line"_Oe6.mpileup > manual_GT/"$line"_Oe6.mpileup.1 ; done
awk '{print $1":"$3, $1,$3,$4,$5}' manual_GT/Step_2.txt | sed 's/ /\t/g' > manual_GT/SNP.bed.2
cat ${SAMPLES} | while read line ; do awk 'NR==FNR{a[$1]=$0 ; next}  {{ found=0; for(i=1;i<=NF;i++) { if($i in a) { print a[$i]; found=1; break; } } if (!found) {print $0} } }' manual_GT/"$line"_Oe6.mpileup.1 manual_GT/SNP.bed.2  > manual_GT/"$line"_Oe6.mpileup.2 ; done

# Step 4: more format modifying 
cat ${SAMPLES} | while read line ; do awk 'BEGIN {OFS="\t"}; { if ($5=="+" ) $5="0"; print $2,$3,$4,$5,$6,$7 }' manual_GT/"$line"_Oe6.mpileup.2 |sed 's/ /\t/g' |   sed 's/\/\A,C,G//g' | sed 's/\/\A,G,C//g' | sed 's/\/\A,C,T//g' | sed 's/\/\A,T,C//g' | sed 's/\/\A,G,T//g' | sed 's/\/\A,T,G//g' | sed 's/\/\C,A,G//g' | sed 's/\/\C,G,A//g' | sed 's/\/\C,A,T//g' | sed 's/\/\C,T,A//g' | sed 's/\/\C,G,T//g' | sed 's/\/\C,T,G//g' | sed 's/\/\G,A,C//g' | sed 's/\/\G,C,A//g' | sed 's/\/\G,A,T//g' | sed 's/\/\G,T,A//g' | sed 's/\/\G,C,T//g' | sed 's/\/\G,T,C//g' | sed 's/\/\T,A,C//g' | sed 's/\/\T,C,A//g' | sed 's/\/\T,A,G//g' | sed 's/\/\T,G,A//g' | sed 's/\/\T,C,G//g' | sed 's/\/\T,G,C//g' | sed 's/\/\A,C//g' | sed 's/\/\C,A//g' | sed 's/\/\A,G//g' | sed 's/\/\G,A//g' | sed 's/\/\A,T//g' | sed 's/\/\T,A//g' | sed 's/\/\C,G//g' | sed 's/\/\G,C//g' | sed 's/\/\C,T//g' | sed 's/\/\T,C//g' | sed 's/\/\G,T//g' | sed 's/\/\T,G//g'  | sed 's/\/\A//g'  | sed 's/\/\C//g' | sed 's/\/\G//g' | sed 's/\/\T//g' > manual_GT/"$line"_Oe6.mpileup.3 ; done

# step 5: again... more format modifying, so that we can run the script countBases2.py to could bases at each position, then delete all mpileup files
cat ${SAMPLES} | while read line ; do awk '{gsub("/","A",$3)}1' manual_GT/"$line"_Oe6.mpileup.3 |  awk '{gsub("+","0",$4)}1' |  awk '{gsub("N","A",$3)}1' | sed 's/ /\t/g'   > manual_GT/"$line"_Oe6.mpileup.4 ; done
cat ${SAMPLES} | while read line ; do python ${countBases2} manual_GT/"$line"_Oe6.mpileup.4 > manual_GT/"$line"_Oe6.count ; done
rm manual_GT/*mpileup*

# step 6: process the count file to add column that contains all the bases called for that position
cat ${SAMPLES} | while read line ; do awk -F "\t" 'BEGIN{OFS="\t"};{print $0, $5+$6+$7+$8+$9+$10}' manual_GT/"$line"_Oe6.count | awk 'BEGIN {OFS="\t"};{max=0;for(i=5;i<9;i++)if($i!~/NA/&&$i>max){max=$i;}; maxfreq=0; if($13==0) maxfreq="-nan" ; else maxfreq=max/$13 ; print $0, maxfreq}' | awk 'BEGIN{OFS="\t"}; { if($14=="-nan" ) $14="0"; print $0}' | awk 'BEGIN{OFS="\t"} {j=sprintf("%8.2f", $14); print $0, j}' | awk 'BEGIN{OFS="\t"}; {a=0; if($13==0) a=0 ; else a=sprintf("%8.2f",$5/$13);print $0, a}' | awk 'BEGIN{OFS="\t"}; {c=0; if($13==0) c=0 ; else c=sprintf("%8.2f",$6/$13);print $0, c}' | awk 'BEGIN{OFS="\t"}; {g=0; if($13==0) g=0 ; else g=sprintf("%8.2f",$7/$13);print $0, g}' | awk 'BEGIN{OFS="\t"}; {t=0; if($13==0) t=0 ; else t=sprintf("%8.2f",$8/$13);print $0, t}' | awk -v col=A 'BEGIN{OFS="\t"}; {NT=""; if($15==$16){NT=col}else{NT="-";}; print $0, NT}' | awk -v col=C 'BEGIN{OFS="\t"}; {NT=""; if($15==$17){NT=col}else{NT="-";}; print $0, NT}' | awk -v col=G 'BEGIN{OFS="\t"}; {NT=""; if($15==$18){NT=col}else{NT="-";}; print $0, NT}' | awk -v col=T 'BEGIN{OFS="\t"}; {NT=""; if($15==$19){NT=col}else{NT="-";}; print $0, NT}'  | awk 'BEGIN{OFS="\t"}; {print $1,$2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $20$21$22$23}' | awk 'BEGIN {OFS="\t"};{if($15=="----")$15="N";print $0}' | awk 'BEGIN {OFS="\t"};{gsub("-","",$15); print $0}' > manual_GT/"$line"_Oe6.count.2 ; done

# step 7: more processing the count file to change the calls based on step 6
cat ${SAMPLES} | while read line ; do awk -v col5=A 'BEGIN {OFS="\t"};{GTA="";if($5!=0){GTA=col5;}else{GTA="-";}print $0,GTA}' manual_GT/"$line"_Oe6.count.2 | awk -v col6=C 'BEGIN {OFS="\t"};{GTC="";if($6!=0){GTC=col6;}else{GTC="-";}print $0,GTC}'| awk -v col7=G 'BEGIN {OFS="\t"};{GTG="";if($7!=0){GTG=col7;}else{GTG="-";}print $0,GTG}' | awk -v col8=T 'BEGIN {OFS="\t"};{GTT="";if($8!=0){GTT=col8;}else{GTT="-";}print $0,GTT}' | awk 'BEGIN {OFS="\t"};{print$0, $16$17$18$19}'| awk 'BEGIN {OFS="\t"};{if($20=="----")$20="N";print $0}'|awk 'BEGIN {OFS="\t"};{gsub("-","",$20); print $0}' | awk 'BEGIN{OFS="\t"}; {print $1,$2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $20}' | awk 'BEGIN {OFS="\t"}; {if($14>0.9 && length($16)>=2) $16=$15; print $0}' | awk 'BEGIN {OFS="\t"}; {if($4>100000000000)$16="N"; print $0}' | awk 'BEGIN {OFS="\t"};{ if(length($16)>=3)$16="N";print$0 }' > manual_GT/"$line"_Oe6.count.3 ; done 

# step 8: extract selected columns into new file and delete previous count files
cat ${SAMPLES} | while read line ; do awk -v sample="$line" 'BEGIN{OFS="\t"}; {print sample, $1":"$2, $16}' manual_GT/"$line"_Oe6.count.3 > manual_GT/"$line"_Oe6.count.4 ; done
rm manual_GT/*count manual_GT/*count.2 manual_GT/*count.3

# step 9: make each SNP diploid (as heterozygous positions will already be)
cat ${SAMPLES} | while read line ; do  awk 'BEGIN{OFS="\t"}; {GT2==""; if($3=="A"){GT2="AA";}else{GT2="-"} print $0, GT2}' manual_GT/"$line"_Oe6.count.4 | awk 'BEGIN{OFS="\t"}; {if($3=="C"){$4="CC";} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="G"){$4="GG";} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="T"){$4="TT";} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="AC"){$4=$3;} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="AG"){$4=$3;} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="AT"){$4=$3;} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="CG"){$4=$3;} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="CT"){$4=$3;} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="GT"){$4=$3;} print $0}' | awk 'BEGIN{OFS="\t"}; {if($3=="N"){$4="NA";} print $1, $2, $4}' > manual_GT/"$line"_Oe6.count.5 ; done

#Step 10: make heterozygous calls, effectively ending up with a vertical fasta file for each position, also remove intermediate files
cat ${SAMPLES} | while read line ; do cat manual_GT/"$line"_Oe6.count.5  | sed 's/"//g' | sed 's/NA/N/g' | sed 's/AC/M/g' | sed 's/AG/R/g' | sed 's/AT/W/g' | sed 's/CG/S/g' | sed 's/CT/Y/g' | sed 's/GT/K/g' | sed 's/AA/A/g' |sed 's/CC/C/g' |sed 's/GG/G/g' | sed 's/TT/T/g' | sed 's/.bam/\n/g' | sed 's/The/>The/g' > manual_GT/"$line"_Oe6_bowtie_unique_f2_sort_manual_GT.vertical.fasta  ; done
rm manual_GT/*count*

# step 11: generate fasta file for each gene, one file per gene, sequences name with sample name
cd manual_GT 
for filename in *.fasta; do tail -n +2 "$filename" | cut -f 3 | tr --delete '\n' | sed 's/^/>'$filename'\n/g' | sed 's/.vertical.fasta//g' | sed '$a\' >> Oe6_bowtie_unique_f2_manual_GT.fa ; done

