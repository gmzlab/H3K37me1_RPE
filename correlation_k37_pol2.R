
############################## LOADING PACKAGES AND FUNCTIONS ##############################

# Setting working directory and loading all necessary packages and data.

setwd('/media/gmzlab/4TB/projects/h3k37/results/')
#install.packages("devtools")
#devtools::install_github("slowkow/ggrepel")
#install.packages("ggpubr")
suppressMessages(library(ggpubr))
#install.packages("ggplot2")
suppressMessages(library(ggplot2))
#BiocManager::install("ggpointdensity")
suppressMessages(library(ggpointdensity))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))


# Loading all annotated features.
#ann <- rtracklayer::import('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens.GRCh38.112.gtf')
# Getting all annotated genes by generating a GRanges object with the corresponding m columns
#write.table(as.matrix(mcols(ann[ann$type == 'gene'])[c('gene_id', 'gene_name')]), '/home/gmzlab/Documents/ricardo/annotation/geneid2name.Homo_sapiens.GRCh38.112.csv', sep = ',', quote = F, row.names = F)
geneid2name <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/geneid2name.Homo_sapiens.GRCh38.112.csv', sep = ',')

quant_groups <- function(x, ngroups = 5, ...) {
  # Calculate quantiles
  quantiles <- quantile(x, probs = seq(0, 1, 1 / ngroups), na.rm = TRUE, ...)
  # Create quantile groups
  groups <- cut(x, breaks = quantiles, include.lowest = TRUE, ...)
  # Print message and return groups
  num.missing <- sum(is.na(groups))
  message(paste('Observations per group: ', paste(table(groups), collapse = ', '),
                '. ', num.missing, ' missing.', sep = ''))
  message('Quantiles: ', levels(groups))  
  levels(groups) <- paste0('Q', 1:ngroups)
  return(as.vector(groups))
}

genes.bed <- read.table("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112_filled.bed")[,1:4]
colnames(genes.bed) <- c("chr", "start", "end", "gene_id", "gene_name", "strand") 
data <- left_join(data, genes.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
promoters.bed <- read.table("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_promoters_1kb_112.bed")[,1:4]
colnames(promoters.bed) <- c("chr", "start", "end", "gene_id") 


############################## LOADING DATA ###############

# Loading tables.
# Log2 over input: marks.
# No normalization: ATAC, RNA Pol II, RNA nascent.

# - Gene bodies using bigwigs: broad marks, ATAC.
data.gene.bodies <- read.table('quantification_deeptools/quantification_gene_bodies.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.gene.bodies)[1] <- "chr"
data.gene.bodies.eu <- read.table('quantification_deeptools/quantification_gene_bodies_eu.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.gene.bodies.eu)[1] <- "chr"
# - Gene bodies using bams: RNA nascent.
data.rna <- read.table('quantification_deeptools/quantification_rna.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.rna)[1] <- "chr"
data.rna.eu <- read.table('quantification_deeptools/quantification_rna_eu.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.rna.eu)[1] <- "chr"
# - Promoters using bigwigs: peak marks, RNA Pol II.
data.promoters <- read.table('quantification_deeptools/quantification_promoters.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.promoters)[1] <- "chr"
data.promoters.eu <- read.table('quantification_deeptools/quantification_promoters_eu.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.promoters.eu)[1] <- "chr"

# Mergin tables into data:
data <- left_join(data.gene.bodies, data.rna, by = c("chr", "start", "end"), relationship = "many-to-many")
data <- left_join(data, genes.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
data.promoters   <- left_join(data.promoters,   promoters.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
data <- left_join(data, data.promoters[,-c(1:3)], by = "gene_id", relationship = "many-to-many")
data <- data[data$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
data <- na.omit(unique(data))

# Mergin tables into data.eu:
data.eu <- left_join(data.gene.bodies.eu, data.rna.eu, by = c("chr", "start", "end"), relationship = "many-to-many")
data.eu <- left_join(data.eu, genes.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
data.promoters.eu   <- left_join(data.promoters.eu,   promoters.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
data.eu <- left_join(data.eu, data.promoters.eu[,-c(1:3)], by = "gene_id", relationship = "many-to-many")
data.eu <- data.eu[data.eu$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
data.eu <- na.omit(unique(data.eu))


rm(data.gene.bodies)
rm(data.promoters)
rm(data.gene.bodies.eu)
rm(data.promoters.eu)
rm(data.rna)
rm(data.rna.eu)


# Operations

data$rna_nascent_1 <- 1e9*data$rna_nascent_1/((data$end - data$start)*as.numeric(system("grep rna_nascent_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$rna_nascent_2 <- 1e9*data$rna_nascent_2/((data$end - data$start)*as.numeric(system("grep rna_nascent_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$rna_nascent <- rowMeans(log2(data[,c("rna_nascent_1", "rna_nascent_2")]+1))
data[,c("rna_nascent_1","rna_nascent_2")] <- NULL
#data[,c(5,12)] <- log2(data[,c(5,12)]+1)
data[,c(4:8,10)] <- log2(data[,c(4:8,10)]+1)

data.eu$rna_nascent_1 <- 1e9*data.eu$rna_nascent_1/((data.eu$end - data.eu$start)*as.numeric(system("grep rna_nascent_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data.eu$rna_nascent_2 <- 1e9*data.eu$rna_nascent_2/((data.eu$end - data.eu$start)*as.numeric(system("grep rna_nascent_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data.eu$rna_nascent <- rowMeans(log2(data.eu[,c("rna_nascent_1", "rna_nascent_2")]+1))
data.eu[,c("rna_nascent_1","rna_nascent_2")] <- NULL
data.eu[,c(4:8,10)] <- log2(data.eu[,c(4:8,10)]+1)


data <- data.gene.bodies.eu

############################## ONE COMPARISON ##############

ptm <- "h3k36me3"

pdf(paste0("plots/scatterplots/scatterplot_", ptm, "_mut_genes_broad_E5.pdf"))
ggplot(data, aes(x = data[, paste0(ptm, "_B12")], y = data[, paste0(ptm, "_B10")])) +
  geom_pointdensity(size = 1, adjust = 2) +
  
  geom_pointdensity(size = 1, adjust = 5, method = "neighbors") +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data[, paste0(ptm, "_B12")]), y = max(data[, paste0(ptm, "_B10")]),
           label = paste0("r = ", round(cor(data[, paste0(ptm, "_B12")], data[, paste0(ptm, "_B10")], method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "WT signal (RPGC)",
         y = "H3.3K37R (RPGC)") +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + coord_cartesian(clip = "off")





dev.off()





pdf("plots/correlations/correlations_pairwise/correlation_eduB12_oris.pdf")
ggplot(data, aes(x = oris, y = edu_B12)) +
  geom_pointdensity(size = 1, adjust = 2) +
  
  geom_pointdensity(size = 1, adjust = 5, method = "kde2d") +
  #geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data$oris), y = max(data$edu_B12),
           label = paste0("r = ", round(cor(data$oris, data$edu_B12, method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "SRR6365075 (log2 CPM)",
         y = "EdU-seq-HU WT (log2 CPM)") +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + coord_cartesian(clip = "off")
dev.off()


############################## SEVERAL COMPARISONS ##############

# Transformar a formato long
data_long <- data %>% pivot_longer(cols = c(oris_1, edu_hal_ctrl_1, edu_hal_nopalb_1), 
                                   names_to = "variable", values_to = "value_x")
data_long$variable <- factor(data_long$variable, levels = c("oris_1", "edu_hal_ctrl_1","edu_hal_nopalb_1"))

# Labels de eje X por faceta
x_labels <- c(oris_1 = "Oris (CPM)",
              edu_hal_ctrl_1 = "Ctrl (CPM)",
              edu_hal_nopalb_1 = "No palb (CPM)")

# Calcular correlación por faceta
r_labels <- data_long %>%
  group_by(variable) %>%
  summarize(
    r = round(cor(value_x, edu_B12, method = "pearson"), 2),
    x_pos = max(value_x, na.rm = TRUE),
    y_pos = max(edu_B12, na.rm = TRUE)
  )

# Crear gráfico
pdf("plots/scatterplots/eduseq.pdf", width= 10, height = 5 )
ggplot(data_long, aes(x = value_x, y = edu_B12)) +
  geom_pointdensity(size = 1, adjust = 20, method = "kde2d") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "#de1f09") +
  geom_text(
    data = r_labels,
    aes(x = x_pos, y = y_pos, label = paste0("r = ", r)),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.1,
    size = 4,
    fontface = "bold", color="#de1f09"
  ) +
  facet_wrap(~variable, scales = "free_x", labeller = labeller(variable = x_labels)) +
  labs(
    y = "B12 (CPM)",
    x = NULL  # el label de eje X ahora lo maneja facet
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12), # título arriba de cada faceta
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  coord_cartesian(clip = "off")
dev.off()


############################## QUANTILES RNA ###########

data$quantile_rna <-as.factor(quant_groups(data$rna_nascent, ngroups = 4))
data$quantile_atac <-as.factor(quant_groups(data$atac, ngroups = 4))
data.eu$quantile_rna <-as.factor(quant_groups(data.eu$rna_nascent, ngroups = 4))

mark <- "h3k37me1_B12"
name <- "H3K37me1"

pdf(paste0("plots/boxplots/quantiles/quantilesrna_",mark,".pdf"))  
ggboxplot(data, x = 'quantile_rna', y = mark, outlier.shape = 16, 
          fill = 'quantile_rna')  + xlab('RNA Quantiles') + #ylim(-2, 2) +
          ylab(paste0(name," (log2 RPGC)")) + theme(legend.position = 'none') +
          scale_fill_brewer(palette='YlOrRd')
dev.off()

pdf(paste0("plots/boxplots/quantiles/quantilesrna_",mark,"_eu.pdf"))  
ggboxplot(data.eu, x = 'quantile_rna', y = mark, outlier.shape = 16, 
          fill = 'quantile_rna')  + xlab('RNA Quantiles') + #ylim(-2, 2) +
          ylab(paste0(name," (log2 RPGC)")) + theme(legend.position = 'none') +
          scale_fill_brewer(palette='YlOrRd')
dev.off()




genes.bed <- read.table("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed")
colnames(genes.bed) <- c("chr", "start", "end", "gene_id","score","strand") 


write.table(left_join(genes.bed,data[c(10,13)], by = "gene_id"), "genes/genes_quantile_rna_nascent.bed", sep = "\t" ,row.names = F, col.names = F, quote = F)



############################## CORRELATION 37 ATAC ###########


data.K37.atac <- read.table('quantification_deeptools/correlation_k37_atac.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.K37.atac)[1] <- "chr"
data.K37.atac$atac <- log2(data.K37.atac$atac+1)
data.K37.atac <- data.K37.atac[data.K37.atac$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
data.K37.atac <- na.omit(unique(data.K37.atac))

pdf("plots/correlations/correlations_pairwise/correlation_pol2_h3k37me1_promoters_gene_bodies.pdf")
ggplot(data.K37.atac, aes(x = atac, y = h3k37me1_async)) +
  geom_pointdensity(size = 1, adjust = 2) +
  
  geom_pointdensity(size = 1, adjust = 2, method = "kde2d") +
  #geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data.K37.atac$atac), y = max(data.K37.atac$h3k37me1_async),
           label = paste0("r = ", round(cor(data.K37.atac$atac, data.K37.atac$h3k37me1_async, method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "ATAC (log2 RPGC)",
         y = "H3K37me1 async (log2 RPGC/input)") +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + coord_cartesian(clip = "off")
dev.off()




############################## LOADING DATA AND FUNCTIONS ##############################

# Setting working directory and functions and loading all necessary packages.
setwd('/media/gmzlab/4TB/projects/h3k37/results/')

suppressMessages(library(tidyverse))
#BiocManager::install("ChIPseeker")
suppressMessages(library(ChIPseeker))
#install.packages("ggupset")
suppressMessages(library(ggupset))
suppressMessages(library(GenomicRanges))
#install.packages("ggpubr")
suppressMessages(library(ggpubr))
#install.packages("cowplot")
suppressMessages(library(cowplot))
#install.packages("heatmaply")
suppressMessages(library(heatmaply))
suppressMessages(library(ggplot2))
#BiocManager::install("ggpointdensity")
suppressMessages(library(ggpointdensity))

# Loading all annotated features.
#ann <- rtracklayer::import('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens.GRCh38.112.gtf')
# Getting all annotated genes by generating a GRanges object with the corresponding m columns
#write.table(as.matrix(mcols(ann[ann$type == 'gene'])[c('gene_id', 'gene_name')]), '/home/gmzlab/Documents/ricardo/annotation/geneid2name.Homo_sapiens.GRCh38.112.csv', sep = ',', quote = F, row.names = F)
geneid2name <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/geneid2name.Homo_sapiens.GRCh38.112.csv', sep = ',')


preTSS <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112_curated_preTSS_5kb_1kb.bed', sep = '\t', header = F)
promoters <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112_curated_promoters_1kb.bed', sep = '\t', header = F)
gene.bodies <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112_curated.bed', sep = '\t', header = F)
gene.bodies.cut <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112_curated_TSS_cut_1kb.bed', sep = '\t', header = F)

colnames(preTSS) <- colnames(promoters) <- colnames(gene.bodies) <- colnames(gene.bodies.cut) <- c("chr", "start","end","gene_id","score","strand")


quant_groups <- function(x, ngroups = 5, ...) {
  # Calculate quantiles
  quantiles <- quantile(x, probs = seq(0, 1, 1 / ngroups), na.rm = TRUE, ...)
  # Create quantile groups
  groups <- cut(x, breaks = quantiles, include.lowest = TRUE, ...)
  # Print message and return groups
  num.missing <- sum(is.na(groups))
  message(paste('Observations per group: ', paste(table(groups), collapse = ', '),
                '. ', num.missing, ' missing.', sep = ''))
  message('Quantiles: ', levels(groups))  
  levels(groups) <- paste0('Q', 1:ngroups)
  return(as.vector(groups))
}
plotting.correlations <- function(xpos, ypos, shift, chip.matrix1, name1, ymin1, ymax1, chip.matrix2, name2, ymin2, ymax2, save=F) {
  
  chip.matrix <- cbind(chip.matrix1[, name1], chip.matrix2[, name2])
  colnames(chip.matrix) <- c(name1, name2)
  
  p1a <- as_tibble(chip.matrix) %>% ggscatter(x = name1, y = name2, add = 'reg.line', size=0.25, 
                                              color='lightgray', conf.int = F, add.params = list(color = 'firebrick', fill = 'lightgray')) +
    stat_cor(label.x = xpos, label.y = ypos) + geom_bin2d(alpha=0.65, bins = 200) +
    scale_fill_continuous(type = 'viridis') + theme(legend.position = 'none') +
    stat_regline_equation(label.x = xpos, label.y = ypos - shift) + xlab(name1) + ylab(name2) +
    xscale('log2', .format = TRUE) + yscale('log2', .format = TRUE) 
  #p1a <- ggpar(p1a, xlim=c(NA, NA), ylim=c(NA, NA)) # ----- DISCARDED
  
  p1b <- as_tibble(chip.matrix) %>% mutate(Quant = as.factor(quant_groups(chip.matrix[,name1]))) %>% 
    mutate(Quant = fct_relevel(Quant, 'Q1','Q2','Q3','Q4','Q5')) %>%
    ggboxplot(x = 'Quant', y = name2, outlier.shape = NA, fill = 'Quant') + ylim(ymin1, ymax1) + 
    xlab(paste0(name1,' Quantiles')) + ylab(name2) + theme(legend.position = 'none') + scale_fill_brewer(palette='YlOrRd')
  
  p1c <- as_tibble(chip.matrix) %>% mutate(Quant = as.factor(quant_groups(chip.matrix[,name2]))) %>% 
    mutate(Quant = fct_relevel(Quant, 'Q1','Q2','Q3','Q4','Q5')) %>%
    ggboxplot(x = 'Quant', y = name1, outlier.shape = NA, fill = 'Quant') + ylim(ymin2, ymax2) +
    xlab(paste0(name2,' Quantiles')) + ylab(name1) + theme(legend.position = 'none') + scale_fill_brewer(palette='YlOrRd')
  
  p2 <-plot_grid(p1a, p1b, p1c, nrow=1, ncol=3, rel_widths=c(2,1.5,1.5), align='hv', axis='tblr')
  
  if (save==T){ ggsave(paste0('plots/correlation_', name1, '_', name2, '.pdf'), width=7, height=5) }
  return(p2)
}

############################## QUANTIFYING COVERAGE SIGNAL ##############################

# Changing file format to SAF:
# - OFS= field separator
system("awk -v OFS='\t' '{ print $4, $1, $2, $3, $6 }' /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_110.bed > ref_all.saf")
system("awk -v OFS='\t' '{ print $4, $1, $2, $3, $6 }' clusters/cluster1_signal.bed > ref_k37.saf")


#intergenic <- read.csv('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_intergenic_110.bed', sep = '\t', header = F)
#intergenic$V4 <- paste0('FP', str_pad(1:nrow(intergenic), 6, 'left', '0'))
#intergenic$V5 <- 0
#intergenic$V6 <- '.'
#write_tsv(intergenic, '/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_intergenic_110.bed', quote = 'none', col_names = F)
system("awk -v OFS='\t' '{ print $4, $1, $2, $3, $6 }' /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_intergenic_110.bed > ref_inter.saf")

# Getting signal counts from .bam files on gene bodies and promoters, specified in .saf files:
# - s: strand-specific (0 = unstranded as default).
# - T: number of threads.
# - F: format of annotation file (.saf).
# - p: paired end reads.
system('featureCounts -s 0 -T 40 -F SAF -p -O -a ref_all.saf -o ip_all_counts.tsv ../bam_files/per_replicate/i*_1.bam')
system('featureCounts -s 0 -T 40 -F SAF -p -O -a ref_all.saf -o all_counts.tsv ../bam_files/per_replicate/edu*.bam ../bam_files/per_replicate/chrna*.bam')
system('featureCounts -s 0 -T 40 -F SAF -p -O -a ref_k37.saf -o k37_counts.tsv ../bam_files/per_replicate/edu*.bam')
system('featureCounts -s 0 -T 40 -F SAF -p -O -a ref_inter.saf -o inter_counts.tsv ../bam_files/per_replicate/edu*.bam')


############################## CHIP-SEQ ##############################
# Importing counts and compute count correlation

tot.counts <- read.csv2('../bam_files/per_replicate/total_counts.csv', header = F, sep = ',', col.names = c('Sample', 'Total'))                
scale.factors <- read.csv2('../bw_files/spikein/scale_factors.txt', header = F, sep = ',', col.names = c('Sample', 'Factor'), row.names = 1, as.is = T)                
scale.factors$Factor <- as.numeric(scale.factors$Factor)

ip.all.counts <- read_tsv('ip_all_counts.tsv', comment = '#') %>% # Reading file without header of comments.
  dplyr::select(-Chr, -Start, -End, -Strand) %>% # Removing some columns.
  gather(-Geneid, -Length, key='Sample', value=Reads) %>% # Reorganizing data.
  mutate(Sample=gsub('../bam_files/per_replicate/','',Sample)) %>% # Removing samples prenames.
  mutate(Sample=gsub('_1.bam','',Sample)) %>% # Removing samples subnames.
  ungroup() %>%
  dplyr::select(Geneid, Sample, Reads) %>%
  group_by(Geneid, Sample) %>% # Sorting data by gene id, strand, length and sample.
  spread(key='Sample', value='Reads')


# Converting gene id column into rownames
ip.all.names <- data.frame(ip.all.counts)[[1]]
ip.all.counts <- data.frame(ip.all.counts[,-1])
row.names(ip.all.counts) <- ip.all.names


# There must not be any value equal to 0.
sum(ip.all.counts[,c('input_B12')]==0)

# Getting spikein-normalized data
ip.all.norm <- ip.all.counts

ip.all.counts['ENSG00000279121',]
ip.all.norm['ENSG00000279121',]
fc.B2['ENSG00000279121']

ip.all.counts['ENSG00000091542',]
ip.all.norm['ENSG00000091542',]
fc.B2['ENSG00000091542']




samples.names <- c('B10', 'B12', 'B2')
for (i in 1:length(samples.names)){
  ip.all.norm[,samples.names[i]] <- ip.all.counts[,paste0('ip_',samples.names[i])]/ip.all.counts[,paste0('input_',samples.names[i])]
  ip.all.norm[,samples.names[i]] <- ip.all.norm[,samples.names[i]]*scale.factors[samples.names[i],'Factor']
}
ip.all.norm <- as.matrix(ip.all.norm[,samples.names])
pdf('/home/gmzlab/Documents/ricardo/lab/conferences/2024_10_workshop_gonzalo/IMAGES_GONZALO/illustrator/boxplot_chipseq_B12_B10.pdf', height = 5 , width =2.5)
boxplot(ip.all.norm[,'B12'],ip.all.norm[,'B10'], outline=F)$stats

dev.off()



boxplot(ip.all.norm[,'B12'],ip.all.norm[,'B10'], outline=F)$stats

wilcox.test(ip.all.norm[,'B12'],ip.all.norm[,'B10'], alternative = 'greater')
median(ip.all.norm[,'B12'])-median(ip.all.norm[,'B10'])

fc.B10 <- log2(ip.all.norm[,'B10']/ip.all.norm[,'B12'])
fc.B2 <- log2(ip.all.norm[,'B2']/ip.all.norm[,'B12'])
boxplot(fc.B10, fc.B2, outline=F)
hist(fc.B10)

positive <- read.csv('clusters/cluster1_signal.bed', header=F,sep='\t') %>% dplyr::select(-V13)

fc.threshold <- -log2(0.7)
#up.B10 <- positive[positive$V4 %in% names(fc.B10[fc.B10 >= fc.threshold]),]
down.B10 <- positive[positive$V4 %in% names(fc.B10[fc.B10 <= -fc.threshold]),]
ns.B10 <- positive[positive$V4 %in% names(fc.B10[fc.B10 > -fc.threshold]),] # & fc.B2 < fc.threshold 
#write_tsv(up.B10, file = 'clusters/up_B10.bed', col_names = F)
write_tsv(down.B10, file = 'clusters/down_B10.bed', col_names = F)
write_tsv(ns.B10, file = 'clusters/ns_B10.bed', col_names = F)

#up.B2 <- positive[positive$V4 %in% names(fc.B2[fc.B2 >= fc.threshold]),]
down.B2 <- positive[positive$V4 %in% names(fc.B2[fc.B2 <= -fc.threshold]),]
ns.B2 <- positive[positive$V4 %in% names(fc.B2[fc.B2 > -fc.threshold]),] # & fc.B2 < fc.threshold 
#write_tsv(up.B2, file = 'clusters/up_B2.bed', col_names = F)
write_tsv(down.B2, file = 'clusters/down_B2.bed', col_names = F)
write_tsv(ns.B2, file = 'clusters/ns_B2.bed', col_names = F)


############################## EDU-SEQ ##############################
tot.counts <- read.csv2('../bam_files/per_replicate/total_counts.tsv', header = T, sep = '\t', col.names = c('Sample', 'Total'))                

k37.counts <- read_tsv('k37_counts.tsv', comment = '#') %>% # Reading file without header of comments.
  dplyr::select(-Chr, -Start, -End, -Strand) %>% # Removing some columns.
  gather(-Geneid, -Length, key='Sample', value=Reads) %>% # Reorganizing data.
  mutate(Sample=gsub('../bam_files/per_replicate/','',Sample)) %>% # Removing samples prenames.
  mutate(Sample=gsub('.bam','',Sample)) %>% # Removing samples subnames.
  left_join(tot.counts, by='Sample') %>%
  mutate(RPKM = (10^9)*Reads/(Total*Length)) %>% 
  ungroup() %>%
  dplyr::select(Geneid, Sample, RPKM) %>%
  mutate(Sample=gsub('edu_','edu',Sample)) %>% # Removing '_' character.
  separate(Sample, into=c('Sample','Rep')) %>% # Getting replicates from samples names.
  group_by(Geneid, Sample) %>% # Sorting data by gene id, strand, length and sample.
  summarise(RPKM = sum(RPKM)/2) %>% # Merging the replicates by count sums.
  spread(key='Sample', value='RPKM')


all.counts <- read_tsv('all_counts.tsv', comment = '#') %>% # Reading file without header of comments.
  dplyr::select(-Chr, -Start, -End, -Strand) %>% # Removing some columns.
  gather(-Geneid, -Length, key='Sample', value=Reads) %>% # Reorganizing data.
  mutate(Sample=gsub('../bam_files/per_replicate/','',Sample)) %>% # Removing samples prenames.
  mutate(Sample=gsub('.bam','',Sample)) %>% # Removing samples subnames.
  left_join(tot.counts, by='Sample') %>%
  mutate(RPKM = (10^9)*Reads/(Total*Length)) %>% 
  ungroup() %>%
  dplyr::select(Geneid, Sample, RPKM) %>%
  mutate(Sample=gsub('edu_','edu',Sample)) %>% # Removing '_' character.
  mutate(Sample=gsub('chrna_','chrna',Sample)) %>% # Removing '_' character.
  separate(Sample, into=c('Sample','Rep')) %>% # Getting replicates from samples names.
  group_by(Geneid, Sample) %>% # Sorting data by gene id, strand, length and sample.
  summarise(RPKM = sum(RPKM)) %>% # Merging the replicates by count sums.
  spread(key='Sample', value='RPKM')

inter.counts <- read_tsv('inter_counts.tsv', comment = '#') %>% # Reading file without header of comments.
  dplyr::select(-Chr, -Start, -End, -Strand) %>% # Removing some columns.
  gather(-Geneid, -Length, key='Sample', value=Reads) %>% # Reorganizing data.
  mutate(Sample=gsub('../bam_files/per_replicate/','',Sample)) %>% # Removing samples prenames.
  mutate(Sample=gsub('.bam','',Sample)) %>% # Removing samples subnames.
  left_join(tot.counts, by='Sample') %>%
  mutate(RPKM = (10^9)*Reads/(Total*Length)) %>% 
  ungroup() %>%
  dplyr::select(Geneid, Sample, RPKM) %>%
  mutate(Sample=gsub('edu_','edu',Sample)) %>% # Removing '_' character.
  separate(Sample, into=c('Sample','Rep')) %>% # Getting replicates from samples names.
  group_by(Geneid, Sample) %>% # Sorting data by gene id, strand, length and sample.
  summarise(RPKM = sum(RPKM)/2) %>% # Merging the replicates by count sums.
  spread(key='Sample', value='RPKM')


# Cleaning up
#system('rm *.tsv.summary *.counts.tsv ref_gene_*')

k37.counts<-as.data.frame(k37.counts)
all.counts<-as.data.frame(all.counts)
inter.counts<-as.data.frame(inter.counts)

rownames(k37.counts) <- k37.counts$Geneid
rownames(all.counts) <- all.counts$Geneid
rownames(inter.counts) <- inter.counts$Geneid

k37.counts <- k37.counts[,-1]
all.counts <- all.counts[,-1]
inter.counts <- inter.counts[,-1]

sum(is.na(inter.counts))

boxplot(all.counts, outline=F)$stats


fc.B10 <- k37.counts[,'eduB10']/k37.counts[,'eduB12']
names(fc.B10) <- rownames(k37.counts)
fc.B2  <- k37.counts[,'eduB2' ]/k37.counts[,'eduB12']
names(fc.B2 ) <- rownames(k37.counts)
boxplot(fc.B10, fc.B2, outline=F)$stats
hist(fc.B10, breaks = 200)
hist(fc.B2 , breaks = 200)

fc.B10 <- all.counts[,'eduB10']/all.counts[,'eduB12']
names(fc.B10) <- rownames(all.counts)
fc.B2  <- all.counts[,'eduB2' ]/all.counts[,'eduB12']
names(fc.B2 ) <- rownames(all.counts)
pdf('plots/fc2.pdf', height = 5 , width =2.5)
boxplot(fc.B2, outline=F, col = 'red3')$stats
dev.off()
hist(fc.B10, breaks = 50)
hist(fc.B2 , breaks = 50)

fc.B10 <- inter.counts[,'eduB10']/inter.counts[,'eduB12']
names(fc.B10) <- rownames(inter.counts)
fc.B2  <- inter.counts[,'eduB2' ]/inter.counts[,'eduB12']
names(fc.B2 ) <- rownames(inter.counts)
boxplot(fc.B10, fc.B2, outline=F)$stats
hist(fc.B10, breaks = 50)
hist(fc.B2 , breaks = 50)





pdf('/home/gmzlab/Documents/ricardo/lab/conferences/2024_10_workshop_gonzalo/IMAGES_GONZALO/illustrator/BOXPLOT_g1_B12_B10.pdf', height = 5 , width =2.5)
boxplot(all.counts[,'eduB12'][fc.B10 <  1],
        all.counts[,'eduB10'][fc.B10 <  1],outline=F, col = c('grey28', 'orangered4'))
dev.off()
wilcox.test(all.counts[,'eduB12'][fc.B10 <  1],
            all.counts[,'eduB10'][fc.B10 <  1], alternative = 'greater')


pdf('/home/gmzlab/Documents/ricardo/lab/conferences/2024_10_workshop_gonzalo/IMAGES_GONZALO/illustrator/BOXPLOT_g2_B12_B10.pdf', height = 5 , width =2.5)


boxplot(all.counts[,'eduB12'][fc.B10 >= 1],
        all.counts[,'eduB10'][fc.B10 >= 1],outline=F, col = c('grey28', 'orangered4'))
dev.off()
wilcox.test(all.counts[,'eduB12'][fc.B10 >=  1],
            all.counts[,'eduB10'][fc.B10 >=  1], alternative = 'less')



pdf('plots/chrna_groups_fc2.pdf', height = 5 , width =2.5)


boxplot(all.counts[,'chrnaB12'][fc.B2 <  1],
        all.counts[,'chrnaB12'][fc.B2 >=  1],outline=F, col=c('#3468a3', '#a3c9ec'))
dev.off()
wilcox.test(all.counts[,'chrnaB12'][fc.B2 <  1],
            all.counts[,'chrnaB12'][fc.B2 >=  1], alternative='two.sided')

pdf('/home/gmzlab/Documents/ricardo/lab/conferences/2024_10_workshop_gonzalo/IMAGES_GONZALO/illustrator/boxplot_chrna_B12_B10.pdf', height = 5 , width =2.5)
boxplot(all.counts[,'chrnaB12'],
        all.counts[,'chrnaB10'],outline=F, col = c('grey28', 'orangered4'))
dev.off()
wilcox.test(all.counts[,'chrnaB12'],
            all.counts[,'chrnaB10'], alternative='two.sided')

pdf('/home/gmzlab/Documents/ricardo/lab/conferences/2024_10_workshop_gonzalo/IMAGES_GONZALO/illustrator/BOXPLOT_EDU_B12.pdf', height = 5 , width =2.5)
boxplot(all.counts[,'eduB12'][fc.B10 <  1],
        all.counts[,'eduB12'][fc.B10 >=  1],outline=F, col = c('grey28', 'orangered4'))
dev.off()
wilcox.test(all.counts[,'eduB12'][fc.B10 <  1],
            all.counts[,'eduB12'][fc.B10 >=  1], alternative='greater')



sum(isNA(log2(k37.counts[,'eduB2' ]/k37.counts[,'eduB12'])))
sum(k37.counts[,'eduB2']==0)
fc.B10['ENSG00000004142']

all.counts[is.na(all.counts[,'eduB2' ]/all.counts[,'eduB12']),]

k37.counts['ENSG00000004897',] == all.counts['ENSG00000004897',]
k37.counts['ENSG00000166226',] == all.counts['ENSG00000166226',]


# ----- LENGTH

# B10
pdf('plots/fc_B10/boxplot_length_downandup.pdf', height = 5 , width =2.5)
boxplot(all.genes[all.genes$name %in% names(fc.B10[fc.B10 <  1])]$width,
        all.genes[all.genes$name %in% names(fc.B10[fc.B10 >= 1])]$width, xlab=c('down', 'up'),outline=F, col=c('#a5272b', '#d19883'))
dev.off()
wilcox.test(all.genes[all.genes$name %in% names(fc.B10[fc.B10 <  1])]$width,
            all.genes[all.genes$name %in% names(fc.B10[fc.B10 >=  1])]$width, alternative = 'less')

pdf('plots/fc_B10/boxplot_length_quartiles1234.pdf', height = 5 , width =2.5)
boxplot(all.genes[all.genes$name %in% names(fc.B10[fc.B10 <  0.9878297])]$width,
        all.genes[all.genes$name %in% names(fc.B10[fc.B10 >= 0.9878297 & fc.B10 <=  1.5235244 ])]$width,
        all.genes[all.genes$name %in% names(fc.B10[fc.B10 >  1.5235244])]$width, xlab=c('q1', 'q2and3', 'q4'),outline=F, col=c('#a5272b', 'white', '#d19883'))
dev.off()

# B2
pdf('plots/fc_B2/boxplot_length_downandup.pdf', height = 5 , width =2.5)
boxplot(all.genes[all.genes$name %in% names(fc.B2[fc.B2 >= 1])]$width,
        all.genes[all.genes$name %in% names(fc.B2[fc.B2 <  1])]$width, xlab=c('up', 'down'),outline=F, col=c('#3468a3', '#a3c9ec'))
dev.off()

pdf('plots/fc_B2/boxplot_length_quartiles1234.pdf', height = 5 , width =2.5)
boxplot(all.genes[all.genes$name %in% names(fc.B2[fc.B2 <  0.8850122])]$width,
        all.genes[all.genes$name %in% names(fc.B2[fc.B2 >= 0.8850122 & fc.B2 <=  1.1138632 ])]$width,
        all.genes[all.genes$name %in% names(fc.B2[fc.B2 >  1.1138632])]$width,xlab=c('q1', 'q2and3', 'q4'),outline=F, col=c('#3468a3', 'white', '#a3c9ec'))
dev.off()

############################## QUARTILES ##############################

fc.B10 <- fc.B10[!is.na(fc.B10)]



names(fc.B10[fc.B10 <  0.9878297])[!is.na(names(fc.B10[fc.B10 <  0.9878297]))]
sum(is.na(fc.B10))




# B10
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B10[fc.B10 <  0.9878297])],                        'genes/fc_B10/fc_quartile1.bed')
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B10[fc.B10 >= 0.9878297 & fc.B10 <= 1.5235244 ])], 'genes/fc_B10/fc_quartile2and3.bed')
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B10[fc.B10 >  1.5235244])],                        'genes/fc_B10/fc_quartile4.bed')
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B10[fc.B10 <  1])], 'genes/fc_B10/fc_quartiledown.bed')
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B10[fc.B10 >= 1])], 'genes/fc_B10/fc_quartileup.bed')




fc.b10.df <-data.frame(name=names(fc.B10), score=fc.B10)

write.table(fc.b10.df[order(fc.b10.df$score, decreasing = T),], 'genes_b10_fc_list.txt', sep = '\t', col.names = T, row.names = F, quote = F)


left_join(as.data.frame(all.genes), fc.b10.df, 'name')



# B2
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B2[fc.B2 <  0.8850122])],                       'genes/fc_B2/fc_quartile1.bed')
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B2[fc.B2 >= 0.8850122 & fc.B2 <= 1.1138632 ])], 'genes/fc_B2/fc_quartile2and3.bed')
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B2[fc.B2 >  1.1138632])],                       'genes/fc_B2/fc_quartile4.bed')
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B2[fc.B2 <  1])], 'genes/fc_B2/fc_quartiledown.bed')
rtracklayer::export(all.genes[all.genes$name %in% names(fc.B2[fc.B2 >= 1])], 'genes/fc_B2/fc_quartileup.bed')


wilcox.test(all.genes[all.genes$name %in% names(fc.B10[fc.B10 > 1.5235244 ])]$width, 
            all.genes[all.genes$name %in% names(fc.B10[fc.B10 < 0.9878297 ])]$width, alternative = 'greater' )



############################## PLOTTING CORRELATIONS ##############################

# Saving new matrices or reloading saved matrices.
save(k37.counts, all.counts, file = '.edu_counts.RData')
load('.edu_counts_.RData') 


plotting.correlations(-4, 4, 1, # X and Y positions of R, p-value and equation text, and Y shift between lines.
                      k37.counts, 'eduB2', 0, 4, # Matrix and name of the first histone ptm and Y limits of the first boxplot.
                      k37.counts, 'eduB12', 0, 4, save=F) # Matrix and name of the second histone ptm and Y limits of the second boxplot.

plotting.correlations(-4, 4, 1, # X and Y positions of R, p-value and equation text, and Y shift between lines.
                      all.counts, 'eduB10', 0, 4, # Matrix and name of the first histone ptm and Y limits of the first boxplot.
                      all.counts, 'eduB2', 0, 4, save=F) # Matrix and name of the second histone ptm and Y limits of the second boxplot.


############################## DISTANCES BETWEEN PEAKS ##############################

setwd('/media/gmzlab/4TB/projects/mut_K37R/results/peaks/')
system('for i in `ls B*.bed |cut -f1 -d _ | sort -u`; do echo $i; bedtools closest -io -iu -D ref -a edu_${i}/${i}_peaks.narrowPeak -b edu_${i}/${i}_peaks.narrowPeak | grep -v "\-1" > ${i}_interpeaks.bed ; done')

B2.interpeaks  <- read.csv2('B2_interpeaks.bed' , sep = '\t', header = F)$V21
B10.interpeaks <- read.csv2('B10_interpeaks.bed', sep = '\t', header = F)$V21
B12.interpeaks <- read.csv2('B12_interpeaks.bed', sep = '\t', header = F)$V21

B2.interpeaks <-B2.interpeaks /max(B2.interpeaks )
B10.interpeaks<-B10.interpeaks/max(B10.interpeaks)
B12.interpeaks<-B12.interpeaks/max(B12.interpeaks)


hist(B2.interpeaks , freq = F, breaks = 200)
hist(B10.interpeaks, breaks = 200)
hist(B12.interpeaks, breaks = 200)
par(mfrow=c(3,1))
hist(B2.interpeaks [B2.interpeaks  < 1000], freq = F, breaks = 200)
hist(B10.interpeaks[B10.interpeaks < 1000], freq = F, breaks = 200)
hist(B12.interpeaks[B12.interpeaks < 1000], freq = F, breaks = 200)
par(mfrow=c(1,1))
boxplot(B2.interpeaks ,
        B10.interpeaks,
        B12.interpeaks, outline=F)

wilcox.test(B10.interpeaks, B12.interpeaks, alternative = 'less')



############################## GETTING GENES WITH A MINIMUM WIDTH ##############################

all.genes <- rtracklayer::import('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_110.bed')

all.genes$width <- width(all.genes)

min.width <- 300000

rtracklayer::export(all.genes[all.genes$width >= min.width,], paste0('genes/all_genes_min_300000.bed'))







############################## VENN DIAGRAM ##############################


vennplot.peakfile(files = c('peaks/edu_B10/B10_peaks.narrowPeak', 'peaks/edu_B12/B12_peaks.narrowPeak'))





############################## MULTIBIGWIGSUMMARY K37 DEPLETION ##############################
sep <- read.table('../results/quantification_deeptools/quantification_h3k37me1_B12_B10.txt',comment.char = "#", header = T)
merge <- sep[,1:3]
merge$H3.3K37R <- rowMeans(sep[,4:5])
merge$WT <- rowMeans(sep[,6:7])



pdf('../results/plots/boxplot_h3k37me1_mut.pdf', height = 5 , width =2.5)
boxplot(merge$WT, merge$H3.3K37R, xlab=c('WT', 'H3.3K37R'),outline=F, col=c('#a5272b', '#d19883'))
dev.off()

wilcox.test(merge$WT, merge$H3.3K37R, alternative = 'two.sided')$p.value


data <- read.delim('quantification_deeptools/b12_vs_halazonetis.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )

data <- read.delim('quantification_deeptools/depletion_fig2_broad_E5.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(data)[1] <- "chr"

data <- na.omit(data[data$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),])
data[,-c(1:3)] <- log2(data[,-c(1:3)] + 1)  
data <- data[rowMeans(data[,4:5]) >0,]

data <- left_join(unique(data), read.csv2("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed", header = F, sep="\t", col.names = c("chr", "start", "end", "gene_id", "score","strand")), by = c("chr", "start", "end"))
data$length <- data$end - data$start

colnames(data) <- gsub("_1","", colnames(data))

data$fc.37  <- data$h3k37me1_B10 / data$h3k37me1_B12
data$fc.h33 <- data$h33_B10 / data$h33_B12

ptm <- "h3k36me3"

pdf(paste0("plots/scatterplots/scatterplot_depletion_", ptm, "_B10vsB12_broad_E5.pdf"))
ggplot(data, aes(x = h3k36me3_B12, y = h3k36me3_B10)) +
  geom_pointdensity(size = 1, adjust = 4, method = "kde2d") +
  geom_abline(slope = 1, intercept = 0, color ="black", linetype="dashed",linewidth = 1) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data[,paste0(ptm,"_B12")]), y = max(data[,paste0(ptm,"_B10")]),
           label = paste0("r = ", round(cor(data[,paste0(ptm,"_B12")], data[,paste0(ptm,"_B10")], method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "WT (log2 RPGC)",
         y = "H3.3K37R (log2 RPGC)") +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + coord_cartesian(clip = "off")
dev.off()

tiff("plots/scatterplots/scatterplot_b12_vs_halazonetis.tiff",
     width = 5,
     height = 5,
     units = "in",
     res = 500,
     compression = "lzw")
pdf(paste0("plots/scatterplots/scatterplot_b12_vs_halazonetis.pdf"), height = 10, width = 10)
ggplot(data, aes(x = oris, y = edu_B12)) +
  geom_pointdensity(size = 1, adjust = 4, method = "kde2d") +
  #geom_abline(intercept = 0, slope = 1, color = "grey20") + 
  #geom_abline(intercept =  0, slope = 1, color = "grey20", linetype="dashed") + 
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = 1.5, y = 0.5,
           label = paste0("r = ", round(cor(data$edu_B12, data$oris_1, method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  #geom_abline(intercept = -log2(fc.threshold), slope = 1, color = "grey20", linetype="dashed") + 
  labs(
    x = "",
    y = ""  ) +
  theme_minimal(base_size = 14) +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + scale_color_gradient(name = "Density") 

dev.off()




############################## MULTIBIGWIGSUMMARY EDU-SEQ ##############################

# -------- LOADING RPGC --------
data <- read.delim('quantification_deeptools/quantification_edu_mut_genes.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
data <- read.delim('quantification_deeptools/quantification_edu_mut_genes_broad_E5.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )

colnames(data)[1] <- "chr"
data <- na.omit(data[data$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),])
data$edu_B10 <- log2(data$edu_B10 + 1)  
data$edu_B12 <- log2(data$edu_B12 + 1)  
data$edu_B2  <- log2(data$edu_B2  + 1)

data$ratio_B10 <- data$edu_B10 - data$edu_B12
data$ratio_B2 <- data$edu_B2 - data$edu_B12
data <- left_join(unique(data), read.csv2("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed", header = F, sep="\t", col.names = c("chr", "start", "end", "gene_id", "score","strand")), by = c("chr", "start", "end"))
data$length <- data$end - data$start

# -------- LOADING COUNTS --------
data <- read.delim('quantification_deeptools/counts_edu_mut_genes_broad_E5.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
data <- read.delim('quantification_deeptools/counts_edu_mut_allgenes.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
data <- read.delim('quantification_deeptools/quantification_edu_mut_rna_nascent_genes_broad_E5_counts.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
data <- read.delim('quantification_deeptools/rna_nascent_genebodies_rawcounts.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(data)[1] <- "chr"
data <- na.omit(data[data$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),])

total.counts.all <- read.delim("/media/gmzlab/4TB/projects/h3k37/bam_files/per_replicate/total_counts.tsv", sep = "\t")
total.counts <- total.counts.all[grep(c("rna_nascent"), total.counts.all$Sample),]
total.counts <- rbind(total.counts,total.counts.all[grep(c("rna.nascent"), total.counts.all$Sample),])




data[,-c(1:3)] <- t(1e6*t(data[,-c(1:3)])/total.counts$Total)
data$rna_nascent_1 <- 1e3*data$rna_nascent_1/(data$end - data$start)
data$rna_nascent_2 <- 1e3*data$rna_nascent_2/(data$end - data$start)
data[,-c(1:3)] <- log2(data[,-c(1:3)]+1)

data$edu_B10 <- rowMeans(cbind(data$edu_B10_1,data$edu_B10_2))
data$edu_B12 <- rowMeans(cbind(data$edu_B12_1,data$edu_B12_2))
data$edu_B2  <- rowMeans(cbind(data$edu_B2_1,data$edu_B2_2))
data$rna_nascent  <- rowMeans(cbind(data$rna_nascent_1,data$rna_nascent_2))


data[,c("edu_B10_1","edu_B10_2","edu_B12_1","edu_B12_2",
        "edu_B2_1","edu_B2_2","rna_nascent_1","rna_nascent_2")] <- NULL


data <- left_join(unique(data), read.csv2("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112_filled.bed", header = F, sep="\t", col.names = c("chr", "start", "end", "gene_id", "score","strand")), by = c("chr", "start", "end"))
data$length <- data$end - data$start
data$ratio_B10 <- data$edu_B10 - data$edu_B12
data$ratio_B2 <- data$edu_B2 - data$edu_B12




# -------- SCATTERPLOT --------

fc.threshold <- 1
data$color_B10 <- ifelse(data$ratio_B10 > log2(fc.threshold), "#E74C3C", ifelse(data$ratio_B10 < -log2(fc.threshold), "#3498DB", "#BDC3C7"))
data$color_B2  <- ifelse(data$ratio_B2  > log2(fc.threshold), "#E74C3C", ifelse(data$ratio_B2  < -log2(fc.threshold), "#3498DB", "#BDC3C7"))

table(data$color_B10)
table(data$color_B2)

pdf(paste0("plots/scatterplots/scatterplot_eduB10_fc", fc.threshold, ".pdf"))
ggplot(data, aes(x = edu_B12, y = edu_B10, color = color_B10)) +
  geom_point(alpha = 1, size = 2) +  # puntos semitransparentes
  #geom_abline(intercept = 0, slope = 1, color = "grey20") + 
  geom_abline(intercept =  log2(fc.threshold), slope = 1, color = "grey20", linetype="dashed") + 
  geom_abline(intercept = -log2(fc.threshold), slope = 1, color = "grey20", linetype="dashed") + 
  geom_smooth(method = "lm", color = "#000000", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data$edu_B12), y = max(data$edu_B10),
           label = paste0("r = ", round(cor(data$edu_B12, data$edu_B10, method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  labs(
    x = "WT log2(RPGC+1)",
    y = "H3.3K37R log2(RPGC+1)"  ) +
  scale_color_identity() +  # Usar colores especificados directamente
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  coord_equal() +  xlim(0, 8.1) +    ylim(0, 8.1)
dev.off()


table(rowSums(data[,c(4,5)] > 0))


pdf(paste0("plots/scatterplots/scatterplot_eduB2_fc", fc.threshold, ".pdf"))
ggplot(data, aes(x = edu_B12, y = edu_B2, color = color_B2)) +
  geom_point(alpha = 1, size = 2) +  # puntos semitransparentes
  geom_abline(intercept = 0, slope = 1, color = "grey20") + 
  geom_abline(intercept =  log2(fc.threshold), slope = 1, color = "grey20", linetype="dashed") + 
  geom_abline(intercept = -log2(fc.threshold), slope = 1, color = "grey20", linetype="dashed") + 
  labs(
    x = "WT log2(RPGC+1)",
    y = "H3.3K37R #2 log2(RPGC+1)"  ) +
  scale_color_identity() +  # Usar colores especificados directamente
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  coord_equal() +  xlim(0, 8.1) +    ylim(0, 8.1)
dev.off()



write.table(data[data$color_B10 =="#E74C3C",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/fc_scatter/groups_E5_B10_up_fc",fc.threshold,".bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
write.table(data[data$color_B10 =="#3498DB",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/fc_scatter/groups_E5_B10_down_fc",fc.threshold,".bed"), quote = F, sep = "\t",row.names = F, col.names = F)
write.table(data[data$color_B2  =="#E74C3C",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/fc_scatter/groups_B2_up_fc",fc.threshold,".bed"),    quote = F, sep = "\t",row.names = F, col.names = F)
write.table(data[data$color_B2  =="#3498DB",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/fc_scatter/groups_B2_down_fc",fc.threshold,".bed"),  quote = F, sep = "\t",row.names = F, col.names = F)


#write.table(data[data$color_B10 =="#E74C3C","gene_id"], paste0("genes/fc_scatter/groups_B10_up_fc",fc.threshold,".txt"),   quote = F, sep = "\t",row.names = F, col.names = F)
#write.table(data[data$color_B10 =="#3498DB","gene_id"], paste0("genes/fc_scatter/groups_B10_down_fc",fc.threshold,".txt"), quote = F, sep = "\t",row.names = F, col.names = F)
#write.table(data[data$color_B2 =="#E74C3C", "gene_id"], paste0("genes/fc_scatter/groups_B2_up_fc",fc.threshold,".txt"),    quote = F, sep = "\t",row.names = F, col.names = F)
#write.table(data[data$color_B2 =="#3498DB", "gene_id"], paste0("genes/fc_scatter/groups_B2_down_fc",fc.threshold,".txt"),  quote = F, sep = "\t",row.names = F, col.names = F)


# -------- GENE LENGTH --------


pdf(paste0('plots/boxplots/groups_edu/gene_length_groups_E4_B10_fc', fc.threshold, '.pdf'), height = 5 , width =2.5)
boxplot(data[data$color_B10 =="#3498DB","length"], data[data$color_B10 =="#E74C3C","length"], 
        names=paste0(c('log2FC < -', 'log2FC > '), fc.threshold),outline=F, col=c('#3498DB', '#E74C3C'))
dev.off()
pdf(paste0('plots/boxplots/groups_edu/gene_length_groups_E4_B2_fc', fc.threshold, '.pdf'), height = 5 , width =2.5)
boxplot(data[data$color_B2 =="#3498DB","length"], data[data$color_B2 =="#E74C3C","length"], 
        names=paste0(c('log2FC < -', 'log2FC > '), fc.threshold),outline=F, col=c('#3498DB', '#E74C3C'))
dev.off()






# -------- QUANTILES --------

data$quant.B12 <- factor(quant_groups(data$edu_B12, ngroups = 3), levels = paste0("Q", 1:3))
data$quant.rna <- factor(quant_groups(data$rna_nascent, ngroups = 4), levels = paste0("Q", 1:4))
#write.table(data[data$quant.B12 =="Q1",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/quartiles_edu_B12_Q1.bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
#write.table(data[data$quant.B12 =="Q2",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/quartiles_edu_B12_Q2.bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
#write.table(data[data$quant.B12 =="Q3",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/quartiles_edu_B12_Q3.bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
#write.table(data[data$quant.B12 =="Q4",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/quartiles_edu_B12_Q4.bed"),   quote = F, sep = "\t",row.names = F, col.names = F)

#   Le cambiamos el nombnre al final
write.table(data[data$quant.rna =="Q1",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/genes_broad_E5_4qs_rna_nascent_Q4.bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
write.table(data[data$quant.rna =="Q2",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/genes_broad_E5_4qs_rna_nascent_Q3.bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
write.table(data[data$quant.rna =="Q3",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/genes_broad_E5_4qs_rna_nascent_Q2.bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
write.table(data[data$quant.rna =="Q4",c("chr", "start", "end", "gene_id", "score","strand")], paste0("genes/genes_broad_E5_4qs_rna_nascent_Q1.bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
table(data$quant.rna)


df.B10 <- data %>%  pivot_longer(cols = c(edu_B12, edu_B10),
                                 names_to = "condition",
                                 values_to = "signal") %>%
  mutate(condition = recode(condition,
                            "edu_B12" = "WT",
                            "edu_B10" = "H3.3K37R"))
df.B10$condition <- factor(df.B10$condition, levels = c("WT", "H3.3K37R"))

pdf(paste0('plots/boxplots/groups_edu/quantiles_rna_B10_cpm.pdf'), height = 5 , width =4)
ggplot(df.B10, aes(x = quant.rna, y = signal, fill = condition)) +
  geom_boxplot(outlier.size = 1, outlier.shape = 16) +
  stat_compare_means(aes(group = condition),
                     method = "wilcox.test",   # O "t.test" si prefieres
                     label = "p.format",
                     label.y = max(df.B10$signal) + 1,  # Ajusta la posición de los p-value
                     size = 2) +  # Ajusta el tamaño del texto
  labs(title = "",
       x = "Transcription quartiles",
       y = "EdU signal (log CPM)") +
  scale_fill_manual(values = c("#a5272b", "#d19883")) +
  theme_bw() + #ylim(0,600)+
  theme(text = element_text(size = 12))
dev.off()

df.B10 %>%
  group_by(quant.B12, condition) %>%
  summarize(mean_signal = mean(signal, na.rm = TRUE)) 


df.B2 <- data %>%  pivot_longer(cols = c(edu_B12, edu_B2),
                                names_to = "condition",
                                values_to = "signal") %>%
  mutate(condition = recode(condition,
                            "edu_B12" = "WT",
                            "edu_B2" = "H3.3K37R #2"))
df.B2$condition <- factor(df.B2$condition, levels = c("WT", "H3.3K37R #2"))

pdf(paste0('plots/boxplots/groups_edu/quantiles_rna_B2_cpm.pdf'), height = 5 , width =4)
ggplot(df.B2, aes(x = quant.rna, y = signal, fill = condition)) +
  geom_boxplot(outlier.size = 1, outlier.shape = 16) +
  stat_compare_means(aes(group = condition),
                     method = "wilcox.test",   # O "t.test" si prefieres
                     label = "p.format",
                     label.y = max(df.B2$signal) + 1,  # Ajusta la posición de los p-value
                     size = 2) +  # Ajusta el tamaño del texto
  labs(title = "",
       x = "Transcription quartiles",
       y = "EdU signal (log RPGC)") +
  scale_fill_manual(values = c("#a5272b", "#6da1c2")) +
  theme_bw() +# ylim(0,9)+
  theme(text = element_text(size = 12))
dev.off()


df.B2 %>%
  group_by(quant.rna, condition) %>%
  summarize(mean_signal = mean(signal, na.rm = TRUE)) 


############################## MULTIBAMSUMMARY ##############################


data <- read.delim('quantification_deeptools/h33_expression.txt', sep = "\t", quote = "'",  header = T )
colnames(data)[1] <- "chr"
data$chr <- as.character(data$chr)
data <- left_join(unique(data), read.csv2("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed", header = F, sep="\t", col.names = c("chr", "start", "end", "gene_id", "score","strand")), by = c("chr", "start", "end"))
data <- left_join(data, geneid2name, "gene_id")

total.counts <- read.delim('../bam_files/per_replicate/total_counts.tsv', sep = "\t", quote = "'",  header = T )
total.counts <-total.counts[total.counts$Sample %in% colnames(data),]

data[,c(4,5)] <- t(1e6*t(data[,c(4,5)])/total.counts$Total)
data[,c(4,5)] <- 1e3*data[,c(4,5)]/(data$end - data$start)
data$rna_nascent <- rowMeans(data[,c(4,5)])

write.table(data[,c("gene_id", "gene_name","rna_nascent_1", "rna_nascent_2", "rna_nascent")], "quantification_deeptools/h33_expression_rpkm.txt",   quote = F, sep = "\t", row.names = F, col.names = T)








############################## MULTIBIGWIGSUMMARY ASYNC VS. MERGE ##############################
sep <- read.delim('quantification_deeptools/pca_async_vs_merge_bs10000.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )

colnames(sep)[1] <- "chr"
data <- sep[,1:3]
data$merge <- rowMeans(sep[,4:5])
data$async <- rowMeans(sep[,6:7])

rm(sep)


data <- na.omit(data[data$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),])
data[,-c(1:3)] <- log2(data[,-c(1:3)] + 1)  

pdf("plots/scatterplots/scatterplot_async_vs_merge_bs10000.pdf")
ggplot(data, aes(x = merge, y = async)) +
  geom_pointdensity(size = 1, adjust = 4, method = "kde2d") +
  geom_abline(slope = 1, intercept = 0, color ="black", linetype="dashed",linewidth = 1) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data[4]), y = max(data[5]),
           label = paste0("r = ", round(cor(data[4], data[5], method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "Merge",
         y = "Asynchronous") +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + coord_cartesian(clip = "off")
dev.off()





############################## MULTIBIGWIGSUMMARY MCM-2 AND RNA ##############################
mcm <- read.delim('quantification_deeptools/mcm2_1_preTSS_5kb_1kb.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
mcm <- read.delim('quantification_deeptools/mcm2_1_promoters_1kb.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
mcm <- read.delim('quantification_deeptools/mcm2_1_TSS_cut_1kb.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(mcm)[1] <- "chr"

mcm <- left_join(mcm, preTSS[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")
mcm <- left_join(mcm, promoters.bed[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")
mcm <- left_join(mcm, gene.bodies.cut[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")


mcm$width <- mcm$end - mcm$start


sum(mcm$width > 29000)
mcm <- mcm[mcm$width > 29000,]

genes <- mcm$gene_id
mcm <- mcm[mcm$gene_id %in% genes,]
#mcm[,4:5] <- mcm[4:5]*1e3 /(mcm$end - mcm$start)



rna <- read.delim('quantification_deeptools/rna_nascent_genebodies_rawcounts.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(rna)[1] <- "chr"
rna$width <- rna$end -rna$start
rna$rna_nascent_1 <- rna$rna_nascent_1 * 1e9 / (rna$width * as.numeric(system("grep rna_nascent_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$rna_nascent_2 <- rna$rna_nascent_2 * 1e9 / (rna$width * as.numeric(system("grep rna_nascent_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$rna <- rowMeans(rna[,4:5])


rna <- read.delim('quantification_deeptools/chrna_B12_genebodies_rawcounts.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
rna <- read.delim('quantification_deeptools/quantification_chrna_B12_B10_allgenes_biomart.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(rna)[1] <- "chr"
rna$width <- rna$end -rna$start
rna$chrna_B12_1 <- rna$chrna_B12_1 * 1e9 / (rna$width * as.numeric(system("grep chrna_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$chrna_B12_2 <- rna$chrna_B12_2 * 1e9 / (rna$width * as.numeric(system("grep chrna_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$rna <- rowMeans(rna[,6:7])




rna <- left_join(rna[,-c(4:8)], genes.bed,  by = c("chr", "start", "end"), relationship = "many-to-many")
rna <- unique(rna)

data <- left_join(mcm[,4:6],rna[,4:5], by="gene_id", relationship = "many-to-many")

data <- data[rowSums(data[,c(1,2,4)]) > 0,]
#data[,c(1,2,4)] <- log2(data[,c(1,2,4)]+1)


data$quartile <- as.factor(quant_groups(data$rna, ngroups = 4))

rna$quartile <- as.factor(quant_groups(rna$rna, ngroups = 4))


table(rna$quartile)
ggplot(data, aes(x = log2(rna+1), y = log2(mcm2_B12_1+1))) +
  geom_pointdensity(size = 1, adjust = 4, method = "kde2d") +
  #geom_abline(slope = 1, intercept = 0, color ="black", linetype="dashed",linewidth = 1) +
  geom_vline(xintercept = log2(0.1+1), color ="#de1f09",linewidth = 0.5) +
  geom_vline(xintercept = log2(0.5+1), color ="#de1f09",linewidth = 0.5) +
  geom_vline(xintercept = log2(2+1), color ="#de1f09",linewidth = 0.5) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  #annotate("text",
  #         x = max(data$rna), y = max(data$mcm2_B12_1),
  #         label = paste0("r = ", round(cor(data$rna, data$mcm2_B12_1, method = "pearson"), 2)),
  #         size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "ChRNA",
         y = "MCM-2") +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + coord_cartesian(clip = "off")



summary(data$rna)
data$quartile <- "no"
data$quartile[data$rna >= 0.1] <- "low"
data$quartile[data$rna >= 0.5] <- "mid"
data$quartile[data$rna >= 1.5] <- "high"
data$quartile <- factor(data$quartile, levels = c("no","low","mid","high"))
table(data$quartile)

# EXPORTAR QUARTILES

quartiles <- na.omit(left_join(gene.bodies.cut, data[,c(3,6)], by = "gene_id"))
write.table(quartiles[quartiles$quartile =="Q1",-7], "genes/genes_4qs_chrna_B12_Q1.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(quartiles[quartiles$quartile =="Q2",-7], "genes/genes_4qs_chrna_B12_Q2.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(quartiles[quartiles$quartile =="Q3",-7], "genes/genes_4qs_chrna_B12_Q3.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(quartiles[quartiles$quartile =="Q4",-7], "genes/genes_4qs_chrna_B12_Q4.bed", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(rna[rna$quartile =="Q1",c(1:3,5:7)], "genes/genes_4qs_chrna_B12_Q4.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(rna[rna$quartile =="Q2",c(1:3,5:7)], "genes/genes_4qs_chrna_B12_Q3.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(rna[rna$quartile =="Q3",c(1:3,5:7)], "genes/genes_4qs_chrna_B12_Q2.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(rna[rna$quartile =="Q4",c(1:3,5:7)], "genes/genes_4qs_chrna_B12_Q1.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# SOLO UNA MUESTRA

summary_data <- data %>%
  group_by(quartile) %>%
  summarise(
    mean_mcm2 = mean(mcm2_B12_1, na.rm = TRUE),
    sem_mcm2 = sd(mcm2_B12_1, na.rm = TRUE) / sqrt(n())
  )

pdf("plots/boxplots/quantiles/mcm2_B12_1_preTSS_5kb_1kb_quartiles_chrna.pdf", height = 5, width = 5)
pdf("plots/boxplots/quantiles/mcm2_B12_1_promoters_1kb_quartiles_chrna.pdf", height = 5, width = 5)
pdf("plots/boxplots/quantiles/mcm2_B12_1_TSS_cut_1kb_quartiles_chrna.pdf", height = 5, width = 5)
ggplot(summary_data, aes(x = quartile, y = mean_mcm2, group = 1)) +
  geom_point(size = 4, color = "#a5272b") +  # Puntos para la media
  geom_line(color = "#a5272b", size = 1, linetype = "dashed") +  # Línea para conectar los puntos
  geom_errorbar(aes(ymin = mean_mcm2 - sem_mcm2, ymax = mean_mcm2 + sem_mcm2), width = 0.2, color = "#a5272b") + # Barras de SEM
  labs(
    x = "Quartile",
    y = "Mean MCM2 Signal (B12_1)",
    title = "Mean MCM2 Signal by Quartile with SEM"
  ) +
  theme_minimal() + ylim(0.8,1.8) +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.5),  # Líneas mayores de la cuadrícula
    axis.line = element_line(color = "black", linewidth = 0.5)  # Líneas de los ejes (X e Y)
  )
dev.off()

shapiro.test(data$mcm2_B12_1[sample(1:length(data$mcm2_B10_1), 5000)])
shapiro.test(data$mcm2_B12_1[sample(1:length(data$mcm2_B12_1), 5000)])

wilcox.test(data$mcm2_B12_1[data$quartile=="Q4"],
            data$mcm2_B12_1[data$quartile=="Q3"])


median(data$mcm2_B10_1[data$quartile=="Q3"]) - median(data$mcm2_B12_1[data$quartile=="Q3"])




# DOS MUESTRAS

long_data <- data %>%
  pivot_longer(cols = c(mcm2_B10_1, mcm2_B12_1), 
               names_to = "sample", 
               values_to = "mcm_signal")

# Calcular media y SEM por cuartil y muestra
summary_data <- long_data %>%
  group_by(quartile, sample) %>%
  summarise(
    mean_signal = mean(mcm_signal, na.rm = TRUE),
    sem_signal = sd(mcm_signal, na.rm = TRUE) / sqrt(n())
  )

# Crear gráfico con puntos y barras de SEM para ambas muestras
pdf("plots/barplots/mcm2_1_preTSS_5kb_1kb_quartiles_rna.pdf", height = 5, width = 5)
pdf("plots/barplots/mcm2_1_promoters_1kb_quartiles_rna.pdf", height = 5, width = 5)
pdf("plots/barplots/mcm2_1_TSS_cut_1kb_quartiles_rna.pdf", height = 5, width = 5)
ggplot(summary_data, aes(x = quartile, y = mean_signal, color = sample, group = sample)) +
  geom_point(size = 4) +  # Puntos para la media
  geom_line(size = 1, linetype = "dashed") +  # Línea para conectar los puntos
  geom_errorbar(aes(ymin = mean_signal - sem_signal, ymax = mean_signal + sem_signal), 
                width = 0.2) +  # Barras de SEM
  labs(
    x = "Quartile",
    y = "Mean MCM2 Signal",
    title = "Mean MCM2 Signal by Quartile for B10 and B12"
  ) +
  scale_color_manual(values = c("#d19883", "#a5272b")) +  # Cambiar colores de las muestras
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_line(color = "gray", size = 0.2),
    axis.line = element_line(color = "black", size = 0.5)
  )
dev.off()

shapiro.test(data$mcm2_B10_1[sample(1:length(data$mcm2_B10_1), 5000)])
shapiro.test(data$mcm2_B12_1[sample(1:length(data$mcm2_B12_1), 5000)])

wilcox.test(data$mcm2_B10_1[data$quartile=="Q2"],
            data$mcm2_B12_1[data$quartile=="Q2"],paired = T)
wilcox.test(data$mcm2_B12_1[data$quartile=="Q4"],
            data$mcm2_B12_1[data$quartile=="Q3"])

median(data$mcm2_B12_1[data$quartile=="Q4"]) - median(data$mcm2_B12_1[data$quartile=="Q3"])
median(data$mcm2_B10_1[data$quartile=="Q3"]) - median(data$mcm2_B12_1[data$quartile=="Q3"])







############################## MULTIBIGWIGSUMMARY MCM-2 B12 VS. B10 ##############################
mcm <- read.delim('quantification_deeptools/mcm2_1_preTSS_5kb_1kb.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(mcm)[c(1,4,5)] <- c("chr","preTSS K37R", "preTSS WT")
mcm <- left_join(mcm, preTSS[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")


mcm.promoters <- read.delim('quantification_deeptools/mcm2_1_promoters_1kb.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(mcm.promoters)[c(1,4,5)] <- c("chr","promoter K37R", "promoter WT")
mcm.promoters <- left_join(mcm.promoters, promoters[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")
mcm <- left_join(mcm[,4:6], mcm.promoters[,4:6],  by = "gene_id", relationship = "many-to-many")
rm(mcm.promoters)

mcm.genebodies <- read.delim('quantification_deeptools/mcm2_1_TSS_cut_1kb.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(mcm.genebodies)[c(1,4,5)] <- c("chr","gene body K37R", "gene body WT")
mcm.genebodies <- left_join(mcm.genebodies, gene.bodies.cut[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")
mcm <- left_join(mcm, mcm.genebodies[,4:6],  by = "gene_id", relationship = "many-to-many")
rm(mcm.genebodies)

mcm <- mcm[,c(3,1,2,4:7)]


mcm <- na.omit(mcm)

mcm <- mcm[rowSums(mcm[,-1]) > 0,]

mcm[,-1] <- log2(mcm[,-1] +1)

mcm$preTSS_FC <- mcm$`preTSS K37R` - mcm$`preTSS WT`
mcm$promoter_FC <- mcm$`promoter K37R` - mcm$`promoter WT`
mcm$gene_body_FC <- mcm$`gene body K37R` - mcm$`gene body WT`


up <- read.csv2('genes/fc_scatter/groups_E5_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_E5_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")

mcm$gene_type <- "ns"
mcm$gene_type[mcm$gene_id %in% up$gene_id] <- "up"
mcm$gene_type[mcm$gene_id %in% down$gene_id] <- "down"







long_data <- mcm %>%
  select(gene_id, preTSS_FC, promoter_FC, gene_body_FC, gene_type) %>%
  pivot_longer(cols = starts_with("preTSS_FC"):starts_with("gene_body_FC"),
               names_to = "region",
               values_to = "FC")

summary_data <- long_data %>%
  group_by(region, gene_type) %>%
  summarise(
    mean_FC = mean(FC, na.rm = TRUE),
    sem_FC = sd(FC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


# Crear el gráfico
ggplot(summary_data, aes(x = region, y = mean_FC, color = gene_type)) +
  geom_point(size = 4) +  # Puntos para la media
  geom_line(aes(group = gene_type), size = 1) +  # Línea para conectar los puntos por tipo de gen
  geom_errorbar(aes(ymin = mean_FC - sem_FC, ymax = mean_FC + sem_FC), width = 0.2) +  # Barras de SEM
  labs(
    x = "Region",
    y = "Fold Change (Mean ± SEM)",
    title = "Fold Change by Region and Gene Type with SEM"
  ) +
  scale_color_manual(values = c("ns" = "gray", "up" = "#E74C3C", "down" = "#3498DB")) +  # Colores para los gene_types
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.5),
    axis.line = element_line(color = "black", size = 0.5)
  )


















############################## MULTIBIGWIGSUMMARY EDU B12 VS. B10 ##############################
edu <- read.delim('quantification_deeptools/edu_preTSS_5kb_1kb_rawcounts.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(edu)[1] <- "chr"
edu$edu_B10_1 <- edu$edu_B10_1 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B10_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$edu_B10_2 <- edu$edu_B10_2 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B10_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$edu_B12_1 <- edu$edu_B12_1 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$edu_B12_2 <- edu$edu_B12_2 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$preTSS_K37R <- rowMeans(edu[,4:5])
edu$preTSS_WT <- rowMeans(edu[,6:7])
edu <- left_join(edu[,-c(4:7)], preTSS[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")


edu.promoters <- read.delim('quantification_deeptools/edu_promoters_1kb_rawcounts.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(edu.promoters)[1] <- "chr"
edu.promoters$edu_B10_1 <- edu.promoters$edu_B10_1 * 1e9 / ((edu.promoters$end -edu.promoters$start) * as.numeric(system("grep edu_B10_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.promoters$edu_B10_2 <- edu.promoters$edu_B10_2 * 1e9 / ((edu.promoters$end -edu.promoters$start) * as.numeric(system("grep edu_B10_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.promoters$edu_B12_1 <- edu.promoters$edu_B12_1 * 1e9 / ((edu.promoters$end -edu.promoters$start) * as.numeric(system("grep edu_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.promoters$edu_B12_2 <- edu.promoters$edu_B12_2 * 1e9 / ((edu.promoters$end -edu.promoters$start) * as.numeric(system("grep edu_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.promoters$promoter_K37R <- rowMeans(edu.promoters[,4:5])
edu.promoters$promoter_WT <- rowMeans(edu.promoters[,6:7])
edu.promoters <- left_join(edu.promoters[,-c(4:7)], promoters[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")

edu <- left_join(edu[,4:6], edu.promoters[,4:6],  by = "gene_id", relationship = "many-to-many")
rm(edu.promoters)

edu.genebodies <- read.delim('quantification_deeptools/edu_TSS_cut_1kb_rawcounts.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(edu.genebodies)[1] <- "chr"
edu.genebodies$edu_B10_1 <- edu.genebodies$edu_B10_1 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B10_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$edu_B10_2 <- edu.genebodies$edu_B10_2 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B10_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$edu_B12_1 <- edu.genebodies$edu_B12_1 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$edu_B12_2 <- edu.genebodies$edu_B12_2 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$gene_body_K37R <- rowMeans(edu.genebodies[,4:5])
edu.genebodies$gene_body_WT <- rowMeans(edu.genebodies[,6:7])
edu.genebodies$width <- edu.genebodies$end -edu.genebodies$start
edu.genebodies <- left_join(edu.genebodies[,-c(4:7)], gene.bodies.cut[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")

edu <- left_join(edu, edu.genebodies[,4:7],  by = "gene_id", relationship = "many-to-many")
rm(edu.genebodies)

#edu <- edu[,c(3,1,2,4:8)]
edu <- edu[,c(3:8)]


edu <- na.omit(edu)

#edu <- edu[rowSums(edu[,2:7]) > 0,]
edu <- edu[rowSums(edu[,2:5]) > 0,]

edu[,2:5] <- log2(edu[,2:5] +1)


edu[,2:7] <- log2(edu[,2:7] +1)

edu$preTSS_FC <- edu$preTSS_K37R - edu$preTSS_WT
edu$promoter_FC <- edu$promoter_K37R - edu$promoter_WT
edu$gene_body_FC <- edu$gene_body_K37R - edu$gene_body_WT

edu$WT_FC <- edu$gene_body_WT - edu$promoter_WT
edu$K37R_FC <- edu$gene_body_K37R - edu$promoter_K37R



up <- read.csv2('genes/fc_scatter/groups_E5_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_E5_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")

edu$gene_type <- "ns"
edu$gene_type[edu$gene_id %in% up$gene_id] <- "up"
edu$gene_type[edu$gene_id %in% down$gene_id] <- "down"




# ----- Grafico de senal -----




long_data <- edu[edu$gene_type =="up" & edu$width > 30000,1:7] %>% gather(-gene_id, key='Sample', value=Signal) %>%
  group_by(gene_id, Sample, Signal)
length(unique(long_data$gene_id))

long_data$genotype <- "WT"
long_data$genotype[long_data$Sample %in% c("preTSS_K37R", "promoter_K37R",  "gene_body_K37R")] <- "H3.3K37R"
long_data$genotype <- factor(long_data$genotype, levels = c("WT", "H3.3K37R"))

long_data$Sample <- gsub("_WT", "", long_data$Sample)
long_data$Sample <- gsub("_K37R", "", long_data$Sample)
long_data$Sample <- factor(long_data$Sample, levels = c("preTSS","promoter",  "gene_body"))


# Crear el gráfico
pdf("plots/boxplots/edu_B12_B10_preTSS_promoters_genebodies_group_up.pdf")
ggplot(long_data, aes(x = Sample, y = Signal, fill = genotype)) +
  geom_boxplot(outliers = FALSE) +  # Quitar los outliers si no los quieres
  stat_compare_means(aes(group = genotype),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 3.5,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) +  # Ajusta el tamaño del texto del p-valor
  labs(title = "EdU-seq-HU Signal Distribution",
       x = "Region",
       y = "EdU Signal (RPKM)") +  
  scale_fill_manual(values = c("#A5272B", "#D19883")) +  # Colores personalizados
  theme_bw() +  # Tema blanco y negro
  theme(text = element_text(size = 12),  # Tamaño de texto general
        axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()


# ----- Grafico de FC entre genotipos ------

edu <- edu[edu$width > 29000,]


df_long <- edu %>%
  filter(gene_type %in% c("up", "down")) %>%
  pivot_longer(
    cols = c(promoter_FC, gene_body_FC),
    names_to = "region",
    values_to = "FC"
  ) %>%
  mutate(
    group = case_when(
      region == "promoter_FC" & gene_type == "up"   ~ "Promoter - up",
      region == "gene_body_FC" & gene_type == "up"  ~ "Gene body - up",
      region == "promoter_FC" & gene_type == "down" ~ "Promoter - down",
      region == "gene_body_FC" & gene_type == "down"~ "Gene body - down"
    ),
    group = factor(group, levels = c(
      "Promoter - up",
      "Gene body - up",
      "Promoter - down",
      "Gene body - down"
    ))
  )

# Colores personalizados

pdf("plots/boxplots/groups_edu/edu_B12_B10_fc_promoters_gene_bodies_min30kb.pdf")
ggplot(df_long, aes(x = group, y = FC, fill = gene_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("up" = "#D55E00",
                               "down" = "#0072B2")) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(    limits = c(-0.9, 0.55),    breaks = seq(-0.75, 0.55, by = 0.25)  )+
  labs(x = "", y = "Log2(H3.3K37R/WT)", fill = "Gene type") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", size = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", size = 0.5),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()


# ----- Grafico de FC entre regiones ------
edu <- edu[edu$width > 29000,]


df_long <- edu %>%
  filter(gene_type %in% c("up", "down")) %>%
  pivot_longer(
    cols = c(WT_FC, K37R_FC),
    names_to = "genotype",
    values_to = "FC"
  ) %>%
  mutate(
    genotype = recode(genotype,
                      WT_FC = "WT",
                      K37R_FC = "K37R"),
    group = case_when(
      genotype == "WT" & gene_type == "up" ~ "WT - UP",
      genotype == "K37R" & gene_type == "up" ~ "K37R - UP",
      genotype == "WT" & gene_type == "down" ~ "WT - DOWN",
      genotype == "K37R" & gene_type == "down" ~ "K37R - DOWN"
    ),
    group = factor(group, levels = c(
      "WT - UP",
      "K37R - UP",
      "WT - DOWN",
      "K37R - DOWN"
    ))
  )



pdf("plots/boxplots/groups_edu/edu_promoters_gene_bodies_fc_B12_B10.pdf")
ggplot(df_long, aes(x = group, y = FC, fill = gene_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("up" = "#D55E00",
                               "down" = "#0072B2")) +
  scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "Gene body / Promoter Fold Change",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", size = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", size = 0.5),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()



df_summary <- df_long %>%
  group_by(group, genotype) %>%
  summarise(
    mean_FC = mean(FC, na.rm = TRUE),
    sem_FC  = sd(FC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


pdf("plots/boxplots/groups_edu/edu_promoters_gene_bodies_fc_B12_B10_sem.pdf")
ggplot(df_summary, aes(x = group, y = mean_FC, color = group)) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(ymin = mean_FC - sem_FC,
        ymax = mean_FC + sem_FC),
    width = 0.2,
    size = 0.8
  ) +
  scale_color_manual(values =c(     "WT - UP" = "#D55E00",
                                    "K37R - UP" = "#f28d3d",
                                    "WT - DOWN"= "#0072B2",
                                    "K37R - DOWN"= "#00a2fc")) +
  scale_y_continuous(    limits = c(-0.8, -0),    breaks = seq(-0.8,0, by = 0.1
                                                               )  )  +
  labs(
    x = "",
    y = "Gene body / Promoter Fold Change",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", size = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", size = 0.5),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()


############################## FIGS 1H Y S1B ##############


data.gene.bodies <- read.table('quantification_deeptools/h3k37me1_atac_gene_bodies_rawcounts.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.gene.bodies)[1] <- "chr"
data.promoters <- read.table('quantification_deeptools/pol2_promoters_rawcounts.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data.promoters)[1] <- "chr"

# Mergin tables into data:
data <- left_join(data.gene.bodies, genes.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
data.promoters   <- left_join(data.promoters,   promoters.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
data.promoters$pol2_1 <- 1e9*data.promoters$pol2_1/((data.promoters$end - data.promoters$start)*as.numeric(system("grep pol2_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data.promoters$pol2_2 <- 1e9*data.promoters$pol2_2/((data.promoters$end - data.promoters$start)*as.numeric(system("grep pol2_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))


data <- left_join(data, data.promoters[,-c(1:3)], by = "gene_id", relationship = "many-to-many")
data <- data[data$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
data <- na.omit(unique(data))


rm(data.gene.bodies)
rm(data.promoters)


# Operations

data$h3k37me1_async_1 <- 1e9*data$h3k37me1_async_1/((data$end - data$start)*as.numeric(system("grep h3k37me1_async_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$h3k37me1_async_2 <- 1e9*data$h3k37me1_async_2/((data$end - data$start)*as.numeric(system("grep h3k37me1_async_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$atac_1 <- 1e9*data$atac_1/((data$end - data$start)*as.numeric(system("grep atac_1 ../bam_files/merge_per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$atac_2 <- 1e9*data$atac_2/((data$end - data$start)*as.numeric(system("grep atac_2 ../bam_files/merge_per_replicate/total_counts.tsv | cut -f2",intern = T)))


data[,c(4:7,9:10)] <- log2(data[,c(4:7,9:10)]+1)

data$h3k37me1_async <- rowMeans(data[,c("h3k37me1_async_1", "h3k37me1_async_2")])
data$atac <- rowMeans(data[,c("atac_1", "atac_2")])
data$pol2 <- rowMeans(data[,c("pol2_1", "pol2_2")])

data <- data[,c(1:3,8,11:13)]







pdf("plots/correlations/correlations_pairwise/correlation_fig1h.pdf")
ggplot(data, aes(x = pol2, y = h3k37me1_async)) +
  #geom_pointdensity(size = 1, adjust = 2) +
  
  geom_pointdensity(size = 1, adjust = 2, method =  "kde2d") +
  #geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data$pol2), y = max(data$h3k37me1_async),
           label = paste0("r = ", round(cor(data$h3k37me1_async, data$pol2, method = "pearson"), 3)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "RNA Pol II ChIP-seq (log2 RPKM)",
         y = "H3K37me1 ChIP-seq (log2 RPKM)", fill = "Density") +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + coord_cartesian(clip = "off")
dev.off()




pdf("plots/correlations/correlations_pairwise/correlation_figs1b.pdf")
ggplot(data, aes(x = atac, y = h3k37me1_async)) +
  #geom_pointdensity(size = 1, adjust = 2) +
  
  geom_pointdensity(size = 1, adjust = 2, method =  "kde2d") +
  #geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data$atac), y = max(data$h3k37me1_async),
           label = paste0("r = ", round(cor(data$h3k37me1_async, data$atac, method = "pearson"), 3)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "ATAC-seq (log2 RPKM)",
         y = "H3K37me1 ChIP-seq (log2 RPKM)", fill = "Density") +
  theme( 
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + coord_cartesian(clip = "off")
dev.off()






############################## FIGS 2C Y S2DE ##############


data <- read.table('quantification_deeptools/quantification_scatters_fig2_rawcount.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data)[1] <- "chr"

data <- left_join(data, genes.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
data <- data[data$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
data <- na.omit(unique(data))



# Operations

for (i in 1:16){
  data[,colnames(data)[4:19][i]]  <- 1e9*data[,colnames(data)[4:19][i]]/((data$end - data$start)*as.numeric(system(paste0("grep ", colnames(data)[4:19][i] , " ../bam_files/per_replicate/total_counts.tsv | grep -v spike | cut -f2"),intern = T)))*as.numeric(system(paste0("grep ", colnames(data)[4:19][i] , " ../bw_files/spikein/scale_factors.txt | cut -f2"),intern = T))
}

data$h33_B10      <-  rowMeans(data[,c("h33_B10_1", "h33_B10_2")])
data$h33_B12      <-  rowMeans(data[,c("h33_B12_1", "h33_B12_2")])
data$h3k36me1_B10 <-  rowMeans(data[,c("h3k36me1_B10_1", "h3k36me1_B10_2")])
data$h3k36me1_B12 <-  rowMeans(data[,c("h3k36me1_B12_1", "h3k36me1_B12_2")])
data$h3k36me3_B10 <-  rowMeans(data[,c("h3k36me3_B10_1", "h3k36me3_B10_2")])
data$h3k36me3_B12 <-  rowMeans(data[,c("h3k36me3_B12_1", "h3k36me3_B12_2")])
data$h3k37me1_B10 <-  rowMeans(data[,c("h3k37me1_B10_1", "h3k37me1_B10_2")])
data$h3k37me1_B12 <-  rowMeans(data[,c("h3k37me1_B12_1", "h3k37me1_B12_2")])





data[,c(4:19)] <- NULL





ptm <- "h33"
ptm <- "h3k36me1"
ptm <- "h3k36me3"
ptm <- "h3k37me1"


############
pdf(paste0("plots/scatterplots/scatterplot_", ptm, "_FINAL.pdf"))
ggplot(data, aes(x = data[, paste0(ptm, "_B12")], y = data[, paste0(ptm, "_B10")])) +
  geom_pointdensity(size = 1, adjust = 5, method = "kde2d") +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", color = "#de1f09", se = T, linewidth = 1, fullrange=T) +
  annotate("text",
           x = 3.9, y = 3.9,
           #x = max(data[, paste0(ptm, "_B12")]), y = max(data[, paste0(ptm, "_B10")]),
           label = paste0("r = ", round(cor(data[, paste0(ptm, "_B12")], data[, paste0(ptm, "_B10")], method = "pearson"), 3)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 17) +
  labs(  x = "Wild type signal (RPKM)", title = ptm,
         y = "H3.3K37R signal (RPKM)") +  
  xlim(c(0,4)) + ylim(c(0,4))     +  
  # ylim(c(0, max(data[, paste0(ptm, c("_B12", "_B10"))]))) +   xlim(c(0, max(data[, paste0(ptm, c("_B12", "_B10"))]))) +
  theme(
    # 🔹 Líneas de los ejes
    axis.line = element_line(color = "black", linewidth = 0.4),
    # 🔹 Líneas guía (ticks)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(4, "pt"),   # largo de los ticks
    # 🔹 Cuadrícula tenue
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) 

dev.off()

############################## FIGS S4BDE ##############

data <- read.table('quantification_deeptools/quantification_mcm2_h3k37me1_genes_figs4abde_allgenes.txt',header = T, sep = "\t",comment.char = "", stringsAsFactors = FALSE)
colnames(data)[1] <- "chr"

data <- left_join(data, genes.bed, by = c("chr", "start", "end"), relationship = "many-to-many")
data <- data[data$chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"),]
data <- unique(data)



# Operations

for (i in 1:8){
  data[,colnames(data)[4:11][i]]  <- 1e9*data[,colnames(data)[4:11][i]]/((data$end - data$start)*as.numeric(system(paste0("grep ", colnames(data)[4:11][i] , " ../bam_files/per_replicate/total_counts.tsv | grep -v spike | cut -f2"),intern = T)))*as.numeric(system(paste0("grep ", colnames(data)[4:11][i] , " ../bw_files/spikein/scale_factors.txt | cut -f2"),intern = T))
}

data$mcm2_B10      <-  rowMeans(data[,c("mcm2_B10_1", "mcm2_B10_2")])
data$mcm2_B12      <-  rowMeans(data[,c("mcm2_B12_1", "mcm2_B12_2")])
data$h3k37me1_B10 <-  rowMeans(data[,c("h3k37me1_B10_1", "h3k37me1_B10_2")])
data$h3k37me1_B12 <-  rowMeans(data[,c("h3k37me1_B12_1", "h3k37me1_B12_2")])





data[,c(4:11)] <- NULL





up <- read.csv2('genes/fc_scatter/groups_allgenes_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_allgenes_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")










data$gene_type <- "ns"
data$gene_type[data$gene_id %in% up$gene_id] <- "Group I"
data$gene_type[data$gene_id %in% down$gene_id] <- "Group II"





wilcox.test(data$h3k37me1_B10[data$gene_type == "Group II"],data$h3k37me1_B12[data$gene_type == "Group II"], paired = F)

wilcox.test(data$mcm2_B10,data$mcm2_B12, paired = T)




df_long <- data %>%
  filter(gene_type %in% c("Group I", "Group II")) %>%
  pivot_longer(
    cols = c(h3k37me1_B12, h3k37me1_B10),
    names_to = "genotype",
    values_to = "Signal"
  ) %>%  mutate(
    genotype = recode(genotype,
                      h3k37me1_B12 = "Wild type",
                      h3k37me1_B10 = "H3.3K37R"))
df_long$genotype <- factor(df_long$genotype, levels = c("Wild type", "H3.3K37R"))



# Crear el gráfico
pdf("plots/boxplots/groups_edu/h3k37me1_B12_B10_groups_allgenes_fc1.pdf", height = 4, width = 5)
ggplot(df_long, aes(x = gene_type, y = Signal, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, staplewidth = 0.5) +
  scale_fill_manual(values = c("Wild type" = "#A5272B",
                               "H3.3K37R" = "#D19883")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "H3K37me1 qChIP-seq signal (RPKM)",
    fill = "Group"
  ) + coord_cartesian(ylim=c(0, 2.2))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
    
  )  +
  stat_compare_means(aes(group = genotype),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 2.1,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) 
dev.off()



df_long <- data %>%
  filter(gene_type %in% c("Group I", "Group II")) %>%
  pivot_longer(
    cols = c(mcm2_B12, mcm2_B10),
    names_to = "genotype",
    values_to = "Signal"
  ) %>%  mutate(
    genotype = recode(genotype,
                      mcm2_B12 = "Wild type",
                      mcm2_B10 = "H3.3K37R"))
df_long$genotype <- factor(df_long$genotype, levels = c("Wild type", "H3.3K37R"))



# Crear el gráfico
pdf("plots/boxplots/groups_edu/mcm2_B12_B10_groups_allgenes_fc1.pdf", height = 4, width = 5)
ggplot(df_long, aes(x = gene_type, y = Signal, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, staplewidth = 0.5) +
  scale_fill_manual(values = c("Wild type" = "#A5272B",
                               "H3.3K37R" = "#D19883")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "MCM2 qChIP-seq signal (RPKM)",
    fill = "Group"
  ) + coord_cartesian(ylim=c(0, 0.5))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
    
  )  +
  stat_compare_means(aes(group = genotype),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 0.49,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) 
dev.off()



df_long <- data  %>%
  pivot_longer(
    cols = c(mcm2_B12, mcm2_B10),
    names_to = "genotype",
    values_to = "Signal"
  ) %>%  mutate(
    genotype = recode(genotype,
                      mcm2_B12 = "Wild type",
                      mcm2_B10 = "H3.3K37R"))
df_long$genotype <- factor(df_long$genotype, levels = c("Wild type", "H3.3K37R"))



# Crear el gráfico
pdf("plots/boxplots/mutants/mcm2_B12_B10_allgenes.pdf", height = 4, width = 4)
ggplot(df_long, aes(x = genotype, y = Signal, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, staplewidth = 0.5) +
  scale_fill_manual(values = c("Wild type" = "#A5272B",
                               "H3.3K37R" = "#D19883")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "MCM2 qChIP-seq signal (RPKM)",
    fill = "Group"
  ) + coord_cartesian(ylim=c(0, 0.5))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
    
  )  +
  stat_compare_means(aes(group = genotype),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 0.49,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) 
dev.off()


