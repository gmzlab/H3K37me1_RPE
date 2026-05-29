#!/bin/bash



############################# LOADING PARAMETERS #############################

# Working directory: absolute path to the folder containing all files.
WD=$(grep working_directory: $1 | awk '{ print $2 }')
# Genome folder absolute path.
GENOME=$(grep genome_path: $1 | awk '{ print $2 }')
# Number of processors.
THREADS=$(grep number_processors: $1 | awk '{ print $2 }')

############################# GETTING GENOMIC BINS COVERAGE #############################



cd $GENOME

#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e 'select chrom, size from hg38.chromInfo' > hg38.pre_genome

#grep -v 'size' hg38.pre_genome > hg38.genome

#rm hg38.pre_genome

#bedtools makewindows -g hg38.genome -w 50000 -s 50000 > genome_50kb.bed

bedtools makewindows -g chromsizes -w 50000 > genome_50kb.bed

cd ${WD}/results

multiBamSummary BED-file -p 40 --BED /home/gmzlab/Documents/ricardo/genome/homo_sapiens_unknown/genome_50kb.bed -b ../bam_files/*.bam -o results.npz --outRawCounts output_rawCount.txt

cat output_rawCount.txt | sort -k1,1 -k2,2n > output_rawCount.sort.txt


############################# CONVERTING BAMS TO BEDS #############################

for i in `ls bam_files_merge/h3k*.bam bam_files_merge/input*.bam | cut -d '/' -f 2 | cut -d '.' -f 1` 
do
	echo $i
	bedtools bamtobed -i bam_files_merge/${i}.bam > bed_files/${i}.bed
done

# When having paired end data, sort by name but do not make any index:
samtools sort -@ 40 -n ${i}.bam -o ${i}_nsorted.bam 
# Then convert into bedpe files:
bedtools bamtobed -i ${i}_nsorted.bam -bedpe > ../../bedpe_files_t2t/${i}.bedpe
# And finally convert into bed files using the start of the forward pair and the end of the reverse end as coordinates of the fragment: 
awk 'BEGIN{OFS="\t"} {print $1, $2, $6, $7, $8, $9}' ${i}.bedpe > ${i}.bed

############################# BINARIZING GENOME COVERAGE #############################

cd ${WD}

# When having a conflict between chromosome names, the next code can be used:
#cd /home/gmzlab/opt/ChromHMM/ANCHORFILES
#for i in `ls * | cut -f1-2 -d .` ; do echo $i ; zcat ${i}.txt.gz | sed 's/^chr//' | gzip > ${i}.fixed.txt.gz ; done
#cd /home/gmzlab/opt/ChromHMM/COORDS
#for i in `ls * | cut -f1-2 -d .` ; do echo $i ; zcat ${i}.bed.gz | sed 's/^chr//' | gzip > ${i}.fixed.bed.gz ; done

java -jar /home/gmzlab/opt/ChromHMM/ChromHMM.jar BinarizeBed -b 50000 /home/gmzlab/Documents/ricardo/genome/homo_sapiens_ucsc/chromsizes bed_files metadata/cell_mark_table_bed.tsv results/binarize

java -jar /home/gmzlab/opt/ChromHMM/ChromHMM.jar BinarizeBam -b 50000 /home/gmzlab/Documents/ricardo/genome/homo_sapiens_ucsc/chromsizes bam_files/merge metadata/cell_mark_table_bam.tsv results/binarize

############################# CHROMATIN STATES #############################

java -jar /home/gmzlab/opt/ChromHMM/ChromHMM.jar LearnModel -p 40 -b 50000 results/binarize results/model 4 hg38

# -noenrich if using a genome assembly without data in ChromHMM scripts folder.

############################# GENES #############################

cd /media/4TB/projects/h3k37/results/chromHMM

for i in `seq 1 5` ; do grep E$i model_bs5000/RPE_5_segments.bed > model_bs5000/RPE_5_segments_E${i}.bed ; done

for i in `seq 1 5` ; do bedtools intersect -u -a /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed -b model_bs5000/RPE_5_segments_E${i}.bed -f 0.5 > ../genes/genes_chromHMM/genes_E${i}.bed ; done


wc -l /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed
# 21571 /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed

wc -l ../genes/genes_chromHMM/genes_E*
#  4107 ../genes/genes_chromHMM/genes_E1.bed
#  1004 ../genes/genes_chromHMM/genes_E2.bed
#  2550 ../genes/genes_chromHMM/genes_E3.bed
# 10929 ../genes/genes_chromHMM/genes_E4.bed
#  1372 ../genes/genes_chromHMM/genes_E5.bed
# 19962 total




for i in `seq 1 5` ; do bedtools intersect -u -a /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112.bed -b RPE_5_segments_E${i}.bed -f 0.5 > ../../genes/genes_chromHMM/genes_broad_E${i}.bed ; done 

wc -l /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112.bed 
# 19130 /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112.bed

wc -l ../../genes/genes_chromHMM/genes_broad_E*
#  3820 ../../genes/genes_chromHMM/genes_broad_E1.bed
#   723 ../../genes/genes_chromHMM/genes_broad_E2.bed
#  1889 ../../genes/genes_chromHMM/genes_broad_E3.bed
#  1295 ../../genes/genes_chromHMM/genes_broad_E4.bed
# 10374 ../../genes/genes_chromHMM/genes_broad_E5.bed
# 18101 total










