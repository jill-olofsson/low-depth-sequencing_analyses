#Step 0: prep
REF=/home/$USER/olive_refgenome_Oe6
raw=/home/$USER/map/

# Step 1: filter the SNPs based on individual coverage median (these are examples)
	#Oe6
cat Oe6_med1-2.txt | while read line ; do bcftools filter -i 'DP>=1 & DP<=2 & QUAL>=20' ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf ; done

cat Oe6_med1-4.txt | while read line ; do bcftools filter -i 'DP>=1 & DP<=4 & QUAL>=20' ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf ; done

cat Oe6_med2-6.txt | while read line ; do bcftools filter -i 'DP>=2 & DP<=6 & QUAL>=20' ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf ; done

cat Oe6_med2-8.txt | while read line ; do bcftools filter -i 'DP>=2 & DP<=8 & QUAL>=20' ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf ; done

cat Oe6_med3-10.txt | while read line ; do bcftools filter -i 'DP>=3 & DP<=10 & QUAL>=20' ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf ; done

cat Oe6_med4-14.txt | while read line ; do bcftools filter -i 'DP>=4 & DP<=14 & QUAL>=20' ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf ; done

cat Oe6_med4-16.txt | while read line ; do bcftools filter -i 'DP>=4 & DP<=16 & QUAL>=20' ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov.vcf.gz > ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf ; done

# Step 2: merge the vcf files
	#compress vcf files and index them
cat samples.txt | while read line ; do bgzip ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf ; done
cat samples.txt | while read line ; do bcftools index ${raw}/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf.gz ; done
	#make the vcf_file.txt
ls ${raw}/2_SNPs/*_Oe6_bowtie_unique_f2_sort_allcov_covfilt.vcf.gz > vcf_files_cov_Oe6.txt
	#merge the files
bcftools merge -m id -l vcf_files_cov_Oe6.txt > merged_Oe6_covfilt.vcf

# Step 3: filter the merged file
vcftools --vcf merged_Oe6_covfilt.vcf --out merged_Oe6_covfilt_allcov_80per_mac3 --max-missing 0.2 --mac 3 --recode
	#make bed file
sed -e 's/chr//' merged_Oe6_covfilt_allcov_80per_mac3.recode.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > merged_Oe6_covfilt_allcov_80per_mac3.bed
