## Step 0: prepare for alignment 
# define location reference
REF=/home/$USER/olive_refgenome_Oe6
# define location to data
data=/home/$USER/cleaned_reads

mkdir map
## Step 1: Map in Bowtie2
cat samples.txt | while read line ; do bowtie2 -x ${REF}/Oe6.scaffolds.fa --no-unal -1  ${data}/"$line"_R1.fastq.gz_filtered_trimmed_trimmed.gz -2 ${data}/"$line"_R2.fastq.gz_filtered_trimmed_trimmed.gz -p 2 -S map/"$line"_Oe6_bowtie.sam ; done

## Step 2: convert sam to bam
cat samples.txt | while read line ; do samtools view -S -b -h map/"$line"_Oe6_bowtie.sam > map/"$line"_Oe6_bowtie.bam ; done

## Step 3: only keeping uniquely aligned reads in proper pairs
mkdir map/1_proper_pairs
cat samples.txt | while read line ; do samtools view -h -bq -f2 map/"$line"_Oe6_bowtie.bam > map/1_proper_pairs/"$line"_Oe6_bowtie_unique_f2.bam ; done

## Step 4: Sort and index bam-file
cat samples.txt | while read line ;  do samtools sort map/1_proper_pairs/"$line"_Oe6_bowtie_unique_f2.bam > map/1_proper_pairs/"$line"_Oe6_bowtie_unique_f2_sort.bam ; done
	# index
cat samples.txt | while read line ; do samtools index map/1_proper_pairs/"$line"_Oe6_bowtie_unique_f2_sort.bam ; done

## Step 5: call SNPs with mpileup and bcftools
mkdir map/2_SNPs
cat samples.txt | while read line ; do samtools mpileup -q 20 -B -uvIf ${REF}/Oe6.scaffolds.fa map/1_proper_pairs/"$line"_Oe6_bowtie_unique_f2_sort.bam | bcftools call -cv -O v - > map/2_SNPs/"$line"_Oe6_bowtie_unique_f2_sort.vcf ; done



