#Step 0: prep
REF=/home/$USER/olive_refgenome_Oe6
raw=/home/$USER/map/

##Step 1: merge and filter
	#compress vcf files and index them
cat samples.txt | while read line ; do bgzip ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort.vcf ; done
cat samples.txt | while read line ; do bcftools index ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort.vcf.gz ; done
	#make the vcf_file.txt
ls map/2_SNPs/*_Oe6_bowtie_unique_f2_sort.vcf.gz > vcf_files_Oe6_raw.txt
	#merge the files
bcftools merge -m id -l vcf_files_Oe6_raw.txt > merged_Oe6_raw.vcf
sed -e 's/chr//'  merged_Oe6_raw.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' >  merged_Oe6_raw.bed

##Step 2: re-call all of the SNPs
cat samples.txt | while read line ; do samtools mpileup -q 20 -B -uvIf ${REF}/Oe6.scaffolds.fa -l merged_Oe6_raw.bed ${raw}/1_proper_pairs/"$line"_Oe6_bowtie_unique_f2_sort.bam | bcftools call -c -O v - > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf ; done

##Step 3: merge and filter
	#compress vcf files and index them
cat samples.txt | while read line ; do bgzip ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf ; done
cat samples.txt | while read line ; do bcftools index ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz ; done
	#make the vcf_file.txt
ls ${raw}/2_SNPs/*_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz > vcf_files_Oe6_allcov.txt
	#merge the files
bcftools merge -m id -l vcf_files_Oe6_allcov.txt > merged_Oe6_allcov.vcf
	#filter the vcf-file
vcftools --vcf merged_Oe6_allcov.vcf --out merged_Oe6_allcov_50per --max-missing  0.5 --recode
	# make a bed file
sed -e 's/chr//' merged_Oe6_allcov_50per.recode.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > merged_Oe6_allcov_50per.bed

##Step 4: re-call all of the SNPs and calc coverages
cat samples.txt | while read line ; do samtools mpileup -q 20 -B -uvIf ${REF}/Oe6.scaffolds.fa -l merged_Oe6_allcov_50per.bed ${raw}/1_proper_pairs/"$line"_Oe6_bowtie_unique_f2_sort.bam | bcftools call -c -O v - > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_50per.vcf ; done
	# get the coverage
cat samples.txt | while read line ; do vcf-query -f '%POS,%QUAL,%INFO/DP\n' ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_50per.vcf | awk -F , '{print $3}' | ./info.sh - >> Oe6_coverag_counts_allcov_50per.txt ; done

