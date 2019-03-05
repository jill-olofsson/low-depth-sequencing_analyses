# Step 0: define location to data
data=/home/$USER/raw_data_tutorial
NGS=/home/$USER/path-to-NGS-QC-Toolkit

# Step 1: NGS QC Toolkit to remove adaptors and remove reads with <80% >q20
mkdir IlluQC
cat samples.txt | while read line ; do perl ${NGS}/QC/IlluQC_PRLL.pl -pe ${data}/"$line"_R1.fastq.gz ${data}/"$line"_R2.fastq.gz 2 A -l 80 -s 20 -c 4 -o IlluQC/ ; done

# Step 3: NGS QC Toolkit to remove and sequence with ambigious base
mkdir AmbiguityFilter 
cat samples.txt | while read line ; do perl ${NGS}/Trimming/AmbiguityFiltering.pl -i IlluQC/"$line"_R1.fastq.gz_filtered -irev IlluQC/"$line"_R2.fastq.gz_filtered -c 0 ; done
mv IlluQC/*_R1.fastq.gz_filtered_trimmed AmbiguityFilter | mv IlluQC/*_R2.fastq.gz_filtered_trimmed AmbiguityFilter  

# Step 4: NGS QC Toolkit to trim low quality bases (<20) from 3' end of sequence
TrimmingReads=/usr/local/extras/Genomics/apps/ngsqcoolkit/current/Trimming
cat samples.txt | while read line ; do perl ${NGS}/Trimming/TrimmingReads.pl -i AmbiguityFilter/"$line"_R1.fastq.gz_filtered_trimmed -irev AmbiguityFilter/"$line"_R2.fastq.gz_filtered_trimmed -q 20 ; done
gzip AmbiguityFilter/*_filtered_trimmed_trimmed

mkdir cleaned_reads
mv AmbiguityFilter/*_R1.fastq.gz_filtered_trimmed_trimmed.gz /cleaned_reads | mv AmbiguityFilter/*_R2.fastq.gz_filtered_trimmed_trimmed.gz /cleaned_reads


