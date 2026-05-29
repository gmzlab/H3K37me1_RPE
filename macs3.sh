#!/bin/bash

cd ${WD}/results/peaks



# Llamada a picos con macs3. La primera vez cutoff-analysis, a partir de ahi ya no.

macs3 callpeak -t ../bam_files/per_replicate/${i}.bam --name $i -f BAMPE --outdir peaks/$i -g hs -q 0.05 --cutoff-analysis




# Seleccion de picos por encima de un valor de fold change y de -log q-value:

cd ${WD}/results/peaks

SAMPLE=mut_vs_parental ; F=4 ; Q=1.3 ; LC_ALL=C sort -k7,7nr ${SAMPLE}/${SAMPLE}_peaks.narrowPeak | LC_ALL=C awk -v t="${F}" '$7 >= t' | LC_ALL=C awk -v t="${Q}" '$9 >= t' > ${SAMPLE}/${SAMPLE}_peaks_fc${F}_q${Q}.narrowPeak ; bedtools intersect -u -a /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112_expanded_1kb.bed -b ${SAMPLE}/${SAMPLE}_peaks_fc${F}_q${Q}.narrowPeak > ${SAMPLE}/${SAMPLE}_genes.bed ; wc -l ${SAMPLE}/${SAMPLE}_peaks_fc${F}_q${Q}.narrowPeak ${SAMPLE}/${SAMPLE}_genes.bed













