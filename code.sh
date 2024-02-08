# Step 1: Download sequencing data using fastq-dump from the SRA toolkit.
# This command splits the files if they were uploaded as paired-end reads.
fastq-dump --split-files SRA2089358

# Step 2: Download the human reference genome (hg38) and decompress it.
# The reference genome is essential for aligning sequencing reads to a known sequence.
wget -P ~/Desktop/bioInf/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ~/Desktop/bioInf/supporting_files/hg38/hg38.fa.gz

# Step 3: Create an index of the reference genome using samtools.
# Indexing allows for efficient access and manipulation of the reference genome.
samtools faidx ~/Desktop/bioInf/supporting_files/hg38/hg38.fa

# Step 4: Generate a .dict file of the reference genome using GATK.
# This dictionary file is used by GATK tools for understanding the structure of the reference genome.
gatk CreateSequenceDictionary R=~/Desktop/bioInf/supporting_files/hg38/hg38.fa O=~/Desktop/bioInf/supporting_files/hg38/hg38.dict

# Step 5: Download known sites files for Base Quality Score Recalibration (BQSR) using wget.
# BQSR improves the accuracy of variant calls by adjusting quality scores of sequencing reads based on known variant sites.
wget -P ~/Desktop/bioInf/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ~/Desktop/bioInf/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# Define paths for input and output directories to organize the workflow.
ref=~/Desktop/bioInf/supporting_files/hg38/hg38.fa
known_sites=~/Desktop/bioInf/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf
aligned_reads=~/Desktop/bioInf/aligned_reads
reads=~/Desktop/bioInf/raw_reads
results=~/Desktop/bioInf/results
data=~/Desktop/bioInf/data

# Quality Control (QC) using FastQC: Assess the quality of raw sequencing reads.
fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# Align reads to reference genome using BWA-MEM: This step aligns sequencing reads to the reference genome to identify genomic coordinates.
bwa index ${ref}
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

# Mark duplicates and sort reads using GATK: Identifies and marks duplicate reads to reduce the impact of PCR duplication artifacts.
gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

# Base Quality Score Recalibration (BQSR) using GATK: Adjusts the quality scores of sequencing reads based on known variant sites and recalibrates base quality scores.
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam

# Collect metrics: This step gathers various metrics to evaluate the quality of the alignment, including alignment summary and insert size metrics.
gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf



# Variant calling using GATK HaplotypeCaller: Calls SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. This step is crucial for identifying genomic variants from aligned reads.

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf

#Variants analysis may be performed using https://www.cancergenomeinterpreter.org
