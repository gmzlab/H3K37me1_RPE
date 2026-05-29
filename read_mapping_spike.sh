#!/bin/bash

# LOADING PARAMETERS:

# Number of processors.
THREADS=$1
# Path to genome file.
GENOME_FILE=$2
# Sample name including input or histone.
SAMPLE=$3
# Wether reads are paired end or not.
PAIRED_END=$4


#echo $THREADS $GENOME_FILE $SAMPLE $PAIRED_END

# MAPPING:

if [ $PAIRED_END == 'yes' ]
then
	
	# Mapping with BWA aligner:
	# -n: maximum edit distance.
	# -k: maximum edit distance in the seed.
	# -R: maximum number of equally best hits before proceeding with suboptimal alignments.
	# -t: number of processors.
	bwa aln -n 3 -k 2 -R 300 -t $THREADS $GENOME_FILE ${SAMPLE}_1.fastq > ${SAMPLE}_spike_1.sai
	bwa aln -n 3 -k 2 -R 300 -t $THREADS $GENOME_FILE ${SAMPLE}_2.fastq > ${SAMPLE}_spike_2.sai
	# -n: maximum number of aligments to output in the XA (alternative hits) tag for reads paired properly.
	bwa sampe -n 3 $GENOME_FILE ${SAMPLE}_spike_1.sai ${SAMPLE}_spike_2.sai ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq > ${SAMPLE}_spike.sam
	
	# Cleaning up:
	# -f: ignore if file does not exist.
	rm -f  ${SAMPLE}_spike_1.sai ${SAMPLE}_spike_2.sai ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq
	  
else
	# Mapping with BWA aligner:
	# -n: maximum edit distance.
	# -k: maximum edit distance in the seed.
	# -R: maximum number of equally best hits before proceeding with suboptimal alignments.
	# -n: maximum number of aligments to output in the XA (alternative hits) tag for reads paired properly.
	bwa aln -n 3 -k 2 -R 300 -t $THREADS $GENOME_FILE ${SAMPLE}.fastq > ${SAMPLE}_spike.sai
	bwa samse -n 3 $GENOME_FILE ${SAMPLE}_spike.sai ${SAMPLE}.fastq > ${SAMPLE}_spike.sam
	
	# Cleaning up:
	# -f: ignore if file does not exist.
	rm -f ${SAMPLE}.fastq ${SAMPLE}_spike.sai

fi

# Filtering unique and high quality mapping reads:
# -F: exclude alignments with the given flag. FLAG 780 refers to unmapped reads, unmapped mates, secondary alignment and reads that fail platform/vendor quality checks (https://broadinstitute.github.io/picard/explain-flags.html).
# -q: exclude reads with mapping quality under a given value. 30 means a 0.001 probability of mismapping (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577856/).
# -O: output format.
# -@: number of processors.
samtools view -F 780 -O BAM -q 30 -@ $THREADS ${SAMPLE}_spike.sam > ${SAMPLE}_spike_smap.bam

# PCR DUPLICATES REMOVAL:

# Sorting .bam file by read name:
# -@: number of processors.
# -n: sort by read names.
# -o: output file.
samtools sort -@ $THREADS -n ${SAMPLE}_spike_smap.bam -o ${SAMPLE}_spike_sorted.bam 

# Filling in mate coordinates and inserting size fields:
# -@: number of processors.
# -m: add mate score tags used later by markdup.
samtools fixmate -@ $THREADS -m ${SAMPLE}_spike_sorted.bam ${SAMPLE}_spike_fixmate.bam

# Sorting again but by coordinates:
# -@: number of processors.
# -o: output file.
samtools sort -@ $THREADS ${SAMPLE}_spike_fixmate.bam -o ${SAMPLE}_spike_rmdup_sorted.bam 

# Removing PCR duplicates:
# -f: file on which save the statistics.
# -@: number of processors.
# -r: remove duplicate reads.
samtools markdup -f ${SAMPLE}_spike_duplicate_remove_stats.txt -@ $THREADS -r ${SAMPLE}_spike_rmdup_sorted.bam ${SAMPLE}_spike.bam

mv ${SAMPLE}_spike_duplicate_remove_stats.txt -t ../../quality_control

echo ''
echo 'Number of reads in original .sam file:'
samtools view -c ${SAMPLE}_spike.sam 
echo 'Number of reads after PCR duplicates removal:'
samtools view -c ${SAMPLE}_spike.bam
echo ''

# Indexing .bam file:
# -@: number of processors.
samtools index -@ $THREADS ${SAMPLE}_spike.bam

# Cleaning up:
# -f: ignore if file does not exist.
rm -f ${SAMPLE}_spike.sam ${SAMPLE}_spike_smap.bam ${SAMPLE}_spike_sorted.bam ${SAMPLE}_spike_rmdup_sorted.bam ${SAMPLE}_spike_fixmate.bam 


