
############################## LOADING DATA AND FUNCTIONS ##############################

# Setting working directory and functions and loading all necessary packages.
setwd('/media/gmzlab/4TB/projects/h3k37/results/')

#install.packages("tidyverse")
suppressMessages(library(tidyverse))
#BiocManager::install("ChIPseeker")
suppressMessages(library(ChIPseeker))
#install.packages("ggupset")
#suppressMessages(library(ggupset))
#suppressMessages(library(GenomicRanges))
#install.packages("ggpubr")
suppressMessages(library(ggpubr))
#install.packages("cowplot")
#suppressMessages(library(cowplot))
#install.packages("heatmaply")
#suppressMessages(library(heatmaply))
suppressMessages(library(ggplot2))
#BiocManager::install("ggpointdensity")
suppressMessages(library(ggpointdensity))
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))

suppressMessages(library(preprocessCore))

suppressMessages(library(ggrepel))



gene.bodies <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112_filled.bed', sep = '\t', header = F)
colnames(gene.bodies) <- c("chr", "start", "end", "gene_id", "gene_name", "strand")

#write.table(gene.bodies[(gene.bodies$end - gene.bodies$start) >= 30000,], "/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112_min30kb.bed", row.names = F, col.names = F, quote = F, sep = "\t")



gene.promoters <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_promoters_1kb_112_min3kb.bed', sep = '\t', header = F)
colnames(gene.promoters) <- c("chr", "start", "end", "gene_id", "gene_name", "strand")

gene.bodies.cut <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_genebodycut_1kb_112_min3kb.bed', sep = '\t', header = F)
colnames(gene.bodies.cut) <- c("chr", "start", "end", "gene_id", "gene_name", "strand")

gene.bodies.3kb <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112_plus3kbTSS.bed', sep = '\t', header = F)
colnames(gene.bodies.3kb) <- c("chr", "start", "end", "gene_id", "gene_name", "strand")

gene.bodies.1kb <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112_plus1kbpreTSS.bed', sep = '\t', header = F)
colnames(gene.bodies.1kb) <- c("chr", "start", "end", "gene_id", "gene_name", "strand")





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


############################## MULTIBAMSUMMARY vs. MULTIBIGWIGSUMMARY ##############################



data <- unique(read.delim('quantification_deeptools/quantification_edu_B12_1_genes_broad_E5.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T ))
colnames(data)[1] <- "chr"
colnames(data)[4] <- "cpm_b12"
colnames(data)[5] <- "cpm_b10"

data.counts <- unique(read.delim('quantification_deeptools/quantification_edu_B12_1_genes_broad_E5_raw.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T ))
colnames(data.counts)[1] <- "chr"
colnames(data.counts)[4] <- "counts_b12"
colnames(data.counts)[5] <- "counts_b10"

data<- left_join(data, data.counts, by = c("chr", "start", "end"))
rm(data.counts)
 
total.counts.b12 <- read.delim("/media/gmzlab/4TB/projects/h3k37/bam_files/per_replicate/total_counts.tsv", sep = "\t")[grep(c("edu_B12_1"), read.delim("/media/gmzlab/4TB/projects/h3k37/bam_files/per_replicate/total_counts.tsv", sep = "\t")$Sample),]$Total
total.counts.b10 <- read.delim("/media/gmzlab/4TB/projects/h3k37/bam_files/per_replicate/total_counts.tsv", sep = "\t")[grep(c("edu_B10_1"), read.delim("/media/gmzlab/4TB/projects/h3k37/bam_files/per_replicate/total_counts.tsv", sep = "\t")$Sample),]$Total

data$counts_b12 <- data$counts_b12*1e6/total.counts.b12
data$counts_b10 <- data$counts_b10*1e6/total.counts.b10






data$length <- data$end - data$start


data$fc.cpm <- data$cpm_b10/data$cpm_b12 
data$fc.counts <- data$counts_b10/data$counts_b12 

data$fc <- data$fc.cpm / data$fc.counts

ggplot(data, aes(x = length, y = fc)) +  geom_point(alpha = 1, size = 2) 
ggplot(data, aes(x = cpm_b12, y = cpm_b10)) +  geom_point(alpha = 1, size = 2) +  geom_abline(intercept =  0, slope = 1, color = "grey20", linetype="dashed") 
ggplot(data, aes(x = counts_b12, y = counts_b10)) +  geom_point(alpha = 1, size = 2)  +  geom_abline(intercept =  0, slope = 1, color = "grey20", linetype="dashed") 



############################## ChRNA ##############################

# Matrix construction
rna <- read.delim('quantification_deeptools/quantification_chrna_B12_B10_allgenes_biomart.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(rna)[1] <- "chr"
rna$width <- rna$end -rna$start
rna$chrna_B12_1 <- rna$chrna_B12_1 * 1e9 / (rna$width * as.numeric(system("grep chrna_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$chrna_B12_2 <- rna$chrna_B12_2 * 1e9 / (rna$width * as.numeric(system("grep chrna_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$chrna_B10_1 <- rna$chrna_B10_1 * 1e9 / (rna$width * as.numeric(system("grep chrna_B10_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$chrna_B10_2 <- rna$chrna_B10_2 * 1e9 / (rna$width * as.numeric(system("grep chrna_B10_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$WT   <- rowMeans(rna[,6:7])
rna$H3.3K37R <- rowMeans(rna[,4:5])
rna <- unique(left_join(rna[,-c(4:7)], gene.bodies,  by = c("chr", "start", "end"), relationship = "many-to-many"))

# Matrix filtering
rna <- rna[rowSums(rna[,5:6]) > 0,]
rna <- rna[rna$width >= 30000,]



ggplot(rna, aes(x = WT, y = H3.3K37R)) +
  geom_pointdensity(size = 1, adjust = 4, method = "kde2d") +
  geom_abline(slope = 1, intercept = 0, color ="black", linetype="dashed",linewidth = 1) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  #annotate("text",
  #         x = max(rna$rna), y = max(rna$mcm2_B12_1),
  #         label = paste0("r = ", round(cor(rna$rna, rna$mcm2_B12_1, method = "pearson"), 2)),
  #         size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "RPKM WT",
         y = "RPKM H3.3K37R") +
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


# Quartiles
rna$quartile <- as.factor(quant_groups(rna$WT, ngroups = 4))

summary(rna$WT)
rna$quartile <- "no"
rna$quartile[rna$WT >= 1] <- "low"
rna$quartile[rna$WT >= 2] <- "mid"
rna$quartile[rna$WT >= 3] <- "high"
#rna$quartile <- factor(rna$quartile, levels = c("no","low","mid","high"))
table(rna$quartile)


# EXPORTAR QUARTILES
quartiles <- na.omit(left_join(gene.bodies.cut, rna[,c(3,6)], by = "gene_id"))
write.table(rna[rna$quartile=="Q1",c(1:3,7,8,9)], "genes/genes_4qs_30kb_chrna_B12_Q1.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(rna[rna$quartile=="Q2",c(1:3,7,8,9)], "genes/genes_4qs_30kb_chrna_B12_Q2.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(rna[rna$quartile=="Q3",c(1:3,7,8,9)], "genes/genes_4qs_30kb_chrna_B12_Q3.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(rna[rna$quartile=="Q4",c(1:3,7,8,9)], "genes/genes_4qs_30kb_chrna_B12_Q4.bed", row.names = F, col.names = F, quote = F, sep = "\t")


# Plotting

long_rna <- rna %>% pivot_longer(cols = c(H3.3K37R, WT), names_to = "Genotype", values_to = "RPKM")
long_rna$Genotype <- factor(long_rna$Genotype, levels = c("WT", "H3.3K37R"))
long_rna$quartile <- factor(long_rna$quartile, levels = c("no", "low", "mid", "high"))


# Crear el gráfico
pdf("plots/boxplots/mutants/chrna_B12_B10_allgenes_biomart_quartiles.pdf", height = 4, width = 6)
ggplot(long_rna, aes(x = quartile, y = RPKM, fill = Genotype)) +
  geom_boxplot(outliers = F) +  # Quitar los outliers si no los quieres
  stat_compare_means(aes(group = Genotype),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 6,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) +  # Ajusta el tamaño del texto del p-valor
  labs(y="ChRNA signal (RPKM)", x = "WT quartiles") +  
  scale_fill_manual(values = c("#A5272B", "#D19883")) +  # Colores personalizados
  theme_bw() +  # Tema blanco y negro
  theme(text = element_text(size = 12),  # Tamaño de texto general
        axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()



up <- read.csv2('genes/fc_scatter/groups_allgenes_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_allgenes_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")

rna$gene_type <- "ns"
rna$gene_type[rna$gene_id %in% up$gene_id] <- "Group I"
rna$gene_type[rna$gene_id %in% down$gene_id] <- "Group II"


sum(rna$WT[rna$gene_type == "Group I"] < 0.1)


df_long <- rna %>%
  filter(gene_type %in% c("Group I","Group II")) %>%
  pivot_longer(
    cols = c(WT),
    names_to = "genotype",
    values_to = "FC"
  ) 

df_long$gene_type <- factor(df_long$gene_type, levels = c("WT", "H3.3K37R"))

pdf("plots/boxplots/groups_edu/chrna_B12_groups_fc1.pdf", width = 4, height = 7)
ggplot(df_long, aes(x = gene_type, y = FC, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("WT" = "#A5272B",
                               "H3.3K37R" = "#D19883")) +
  stat_compare_means(aes(group = gene_type),
                     method = "wilcox.test",   # O "t.test" si prefieres
                     label = "p.format",
                     label.y = 4,  # Ajusta la posición de los p-value
                     size = 4) +  # Ajusta el tamaño del texto
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "ChRNA-seq",
    fill = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + coord_cartesian(ylim = c(0, 5))
dev.off()

############################## EdU PEAKS IN REPLI-SEQ REGIONS ##############################

system("grep Early clusters/repli_clusters.bed > clusters/repli_clusters_Early.bed")
system("grep Mid clusters/repli_clusters.bed > clusters/repli_clusters_Mid.bed")
system("grep Late clusters/repli_clusters.bed > clusters/repli_clusters_Late.bed")




df <- data.frame(Time=c("Early", "Mid", "Late"))
df$Constitutive<- c(as.numeric(system("bedtools intersect -u -a peaks/chrom_curated/edu_B12_peaks.bed -b clusters/repli_clusters_Early.bed | wc -l", intern = T)),
                as.numeric(system("bedtools intersect -u -a peaks/chrom_curated/edu_B12_peaks.bed -b clusters/repli_clusters_Mid.bed  | wc -l", intern = T)),
                as.numeric(system("bedtools intersect -u -a peaks/chrom_curated/edu_B12_peaks.bed -b clusters/repli_clusters_Late.bed  | wc -l", intern = T)))

df$K37R_dependent <- c(as.numeric(system("bedtools intersect -u -a peaks/chrom_curated/edu_B10_noB12_peaks.bed -b clusters/repli_clusters_Early.bed  | wc -l", intern = T)),
as.numeric(system("bedtools intersect -u -a peaks/chrom_curated/edu_B10_noB12_peaks.bed -b clusters/repli_clusters_Mid.bed  | wc -l", intern = T)),
as.numeric(system("bedtools intersect -u -a peaks/chrom_curated/edu_B10_noB12_peaks.bed -b clusters/repli_clusters_Late.bed  | wc -l", intern = T)))

df[2:3] <- t(100*t(df[2:3])/colSums(df[2:3]))

pdf("plots/piecharts/edu_constitutive_peaks_in_repli_regions.pdf")
ggplot(df, aes(x = "", y = Constitutive, fill = Time)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  geom_text(aes(label = paste0(Time, "\n", round(Constitutive,2), "%")), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("#429A48",  "#e39610",  "#4f95ff" )) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

pdf("plots/piecharts/edu_K37Rdependent_peaks_in_repli_regions.pdf")
ggplot(df, aes(x = "", y = K37R_dependent, fill = Time)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  geom_text(aes(label = paste0(Time, "\n", round(K37R_dependent, 2), "%")), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("#429A48",  "#e39610",  "#4f95ff" )) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()






############################## EdU GROUPS ##############################

data <- read.delim('quantification_deeptools/quantification_edu_mut_genes_broad_E5.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(data)[1] <- "chr"

data <- read.delim('quantification_deeptools/quantification_edu_mut_allgenes.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(data)[1] <- "chr"


total.counts.all <- read.delim("/media/gmzlab/4TB/projects/h3k37/bam_files/per_replicate/total_counts.tsv", sep = "\t")
total.counts <- total.counts.all[grep(c("edu"), total.counts.all$Sample),]

data[,-c(1:3)] <- t(1e6*t(data[,-c(1:3)])/total.counts$Total)
data[,-c(1:3)] <- 1e3*data[,-c(1:3)]/(data$end - data$start)


data$edu_B10 <- rowMeans(cbind(data$edu_B10_1,data$edu_B10_2))
data$edu_B12 <- rowMeans(cbind(data$edu_B12_1,data$edu_B12_2))
data$edu_B2  <- rowMeans(cbind(data$edu_B2_1,data$edu_B2_2))

data[,c("edu_B10_1","edu_B10_2","edu_B12_1","edu_B12_2",
        "edu_B2_1","edu_B2_2")] <- NULL


data <- left_join(unique(data), gene.bodies, by = c("chr", "start", "end"))
data$length <- data$end - data$start
data$ratio_B10 <- log2(data$edu_B10+1) - log2(data$edu_B12+1)
data$ratio_B2 <- log2(data$edu_B2+1) - log2(data$edu_B12+1)






fc.threshold <- 1.2
data$color_B10 <- ifelse(data$ratio_B10 > log2(fc.threshold), "#E74C3C", ifelse(data$ratio_B10 < -log2(fc.threshold), "#3498DB", "#BDC3C7"))
data$color_B2  <- ifelse(data$ratio_B2  > log2(fc.threshold), "#E74C3C", ifelse(data$ratio_B2  < -log2(fc.threshold), "#3498DB", "#BDC3C7"))

table(data$color_B10)
table(data$color_B2)

highlight <- c("ENSG00000000419", "ENSG00000173338", "ENSG00000142186")

pdf(paste0("plots/scatterplots/scatterplot_eduB10_fc", fc.threshold, ".pdf"))
ggplot(data, aes(x = log2(edu_B12+1), y = log2(edu_B10+1), color = color_B10)) +
  geom_point(alpha = 1, size = 2) +  # puntos semitransparentes
  geom_point(data = subset(data, gene_id %in% highlight), color = "black",size = 2) +
  geom_text_repel(data = subset(data, gene_id %in% highlight),
                  aes(label = gene_name),
                  size = 4, max.overlaps = 10, colour = 'black') +
  #geom_abline(intercept = 0, slope = 1, color = "grey20") + 
  geom_abline(intercept =  log2(fc.threshold), slope = 1, color = "grey20", linetype="dashed") + 
  geom_abline(intercept = -log2(fc.threshold), slope = 1, color = "grey20", linetype="dashed") + 
  geom_smooth(method = "lm", color = "#000000", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(log2(data$edu_B12+1)), y = max(log2(data$edu_B10+1)),
           label = paste0("r = ", round(cor(log2(data$edu_B12+1), log2(data$edu_B10+1), method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  labs(
    x = "WT log2(RPKM+1)",
    y = "H3.3K37R log2(RPKM+1)"  ) +
  scale_color_identity() +  # Usar colores especificados directamente
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank()
  ) +
  coord_equal() #+  xlim(0, 8.1) +    ylim(0, 8.1)
dev.off()




write.table(data[data$color_B10 =="#E74C3C",c("chr", "start", "end", "gene_id", "gene_name","strand")], paste0("genes/fc_scatter/groups_E5_B10_up_fc",fc.threshold,".bed"),   quote = F, sep = "\t",row.names = F, col.names = F)
write.table(data[data$color_B10 =="#3498DB",c("chr", "start", "end", "gene_id", "gene_name","strand")], paste0("genes/fc_scatter/groups_E5_B10_down_fc",fc.threshold,".bed"), quote = F, sep = "\t",row.names = F, col.names = F)


############################## EdU FC GENE BODY - PROMOTER ##############################

edu <- read.delim('quantification_deeptools/quantification_edu_mut_promoters_1kb_min3kb.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(edu)[1] <- "chr"
edu$edu_B10_1 <- edu$edu_B10_1 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B10_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$edu_B10_2 <- edu$edu_B10_2 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B10_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$edu_B12_1 <- edu$edu_B12_1 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$edu_B12_2 <- edu$edu_B12_2 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$edu_B2_1 <- edu$edu_B2_1 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B2_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$edu_B2_2 <- edu$edu_B2_2 * 1e9 / ((edu$end -edu$start) * as.numeric(system("grep edu_B2_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu$promoter_H3.3K37R <- rowMeans(edu[,4:5])
edu$promoter_WT <- rowMeans(edu[,6:7])
edu$promoter_H3.3K37R2 <- rowMeans(edu[,8:9])

edu <- left_join(edu[,-c(4:9)], gene.promoters[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")



edu.genebodies <- read.delim('quantification_deeptools/quantification_edu_mut_genebodycut_1kb_min3kb.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(edu.genebodies)[1] <- "chr"
edu.genebodies$edu_B10_1 <- edu.genebodies$edu_B10_1 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B10_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$edu_B10_2 <- edu.genebodies$edu_B10_2 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B10_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$edu_B12_1 <- edu.genebodies$edu_B12_1 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$edu_B12_2 <- edu.genebodies$edu_B12_2 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$edu_B2_1 <- edu.genebodies$edu_B2_1 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B2_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$edu_B2_2 <- edu.genebodies$edu_B2_2 * 1e9 / ((edu.genebodies$end -edu.genebodies$start) * as.numeric(system("grep edu_B2_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
edu.genebodies$gene_body_H3.3K37R <- rowMeans(edu.genebodies[,4:5])
edu.genebodies$gene_body_WT <- rowMeans(edu.genebodies[,6:7])
edu.genebodies$gene_body_H3.3K37R2 <- rowMeans(edu.genebodies[,8:9])
edu.genebodies$width <- edu.genebodies$end -edu.genebodies$start
edu.genebodies <- left_join(edu.genebodies[,-c(4:9)], gene.bodies.cut[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")

edu <- left_join(edu, edu.genebodies[,4:8],  by = "gene_id", relationship = "many-to-many")
rm(edu.genebodies)

edu <- edu[,c(7,11,4:6,8:10)]
edu <- unique(edu)



#edu <- edu[rowSums(edu[,2:7]) > 0,]

sum(rowSums(edu[,3:8] > 0) <1)

edu <- edu[rowSums(edu[,3:8] > 0) >0,]




edu[,3:8] <- log2(edu[,3:8] +1)




edu$WT_FC <- edu$gene_body_WT / edu$promoter_WT
edu$H3.3K37R_FC <- edu$gene_body_H3.3K37R / edu$promoter_H3.3K37R
edu$H3.3K37R2_FC <- edu$gene_body_H3.3K37R2 / edu$promoter_H3.3K37R2



up <- read.csv2('genes/fc_scatter/groups_allgenes_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_allgenes_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")

edu$gene_type <- "ns"
edu$gene_type[edu$gene_id %in% up$gene_id] <- "up"
edu$gene_type[edu$gene_id %in% down$gene_id] <- "down"




edu <- edu[edu$width > 29000,]

##### B10

median(edu$H3.3K37R_FC[edu$gene_type == "up"])
median(edu$WT_FC[edu$gene_type == "up"])

wilcox.test(edu$H3.3K37R_FC[edu$gene_type == "down"],edu$WT_FC[edu$gene_type == "down"], paired = T)

df_long <- edu %>%
  filter(gene_type %in% c("down", "up")) %>%
  pivot_longer(
    cols = c(WT_FC, H3.3K37R_FC),
    names_to = "genotype",
    values_to = "FC"
  ) %>%
  mutate(
    genotype = recode(genotype,
                      WT_FC = "WT",
                      H3.3K37R_FC = "H3.3K37R"),
    gene_type = recode(gene_type, up="Group I", down = "Group II")
  )

df_long$genotype <- factor(df_long$genotype, levels = c("WT", "H3.3K37R"))

pdf("plots/boxplots/groups_edu/edu_promoters_gene_bodies_fc_B12_B10.pdf", width = 7, height = 5)
ggplot(df_long, aes(x = gene_type, y = FC, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, staplewidth = 0.5) +
  scale_fill_manual(values = c("WT" = "#A5272B",
                               "H3.3K37R" = "#D19883")) +
  stat_compare_means(aes(group = genotype),
                     method = "wilcox.test",   # O "t.test" si prefieres
                     label = "p.format",
                     label.y = 2.1,  # Ajusta la posición de los p-value
                     size = 4) +  # Ajusta el tamaño del texto
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "Fold Change (Gene body / Promoter)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + coord_cartesian(ylim = c(0, 2.2))
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
  scale_y_continuous(    limits = c(-0.9, -0.1),    breaks = seq(-0.8,0, by = 0.1  )  )  +
  labs(
    x = "",
    y = "Log2FC Gene body / Promoter",
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

#### B2

df_long <- edu %>%
  filter(gene_type %in% c("up", "down")) %>%
  pivot_longer(
    cols = c(WT_FC, K37R2_FC),
    names_to = "genotype",
    values_to = "FC"
  ) %>%
  mutate(
    genotype = recode(genotype,
                      WT_FC = "WT",
                      K37R2_FC = "K37R2"),
    group = case_when(
      genotype == "WT" & gene_type == "up" ~ "WT - UP",
      genotype == "K37R2" & gene_type == "up" ~ "K37R2 - UP",
      genotype == "WT" & gene_type == "down" ~ "WT - DOWN",
      genotype == "K37R2" & gene_type == "down" ~ "K37R2 - DOWN"
    ),
    group = factor(group, levels = c(
      "WT - UP",
      "K37R2 - UP",
      "WT - DOWN",
      "K37R2 - DOWN"
    ))
  )



pdf("plots/boxplots/groups_edu/edu_promoters_gene_bodies_fc_B12_B2.pdf")
ggplot(df_long, aes(x = group, y = FC, fill = gene_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("up" = "#D55E00",
                               "down" = "#0072B2")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "Log2FC Gene body / Promoter",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.5),
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


pdf("plots/boxplots/groups_edu/edu_promoters_gene_bodies_fc_B12_B2_sem.pdf")
ggplot(df_summary, aes(x = group, y = mean_FC, color = group)) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(ymin = mean_FC - sem_FC,
        ymax = mean_FC + sem_FC),
    width = 0.2,
    size = 0.8
  ) +
  scale_color_manual(values =c(     "WT - UP" = "#D55E00",
                                    "K37R2 - UP" = "#f28d3d",
                                    "WT - DOWN"= "#0072B2",
                                    "K37R2 - DOWN"= "#00a2fc")) +
  scale_y_continuous(    limits = c(-0.9, -0.1),    breaks = seq(-0.8,0, by = 0.1  )  )  +
  labs(
    x = "",
    y = "Log2FC Gene body / Promoter",
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






############################## EdU BOXPLOTS ##############################

#table <- read.csv2("genes/genes_chromHMM/genes_broad_E5.bed", sep = "\t", header = F)
#table[table$V5 == "","V5"] <- "NA"
#write.table(table, "genes/genes_chromHMM/genes_broad_E5_filled.bed", quote = F, sep = "\t",row.names = F, col.names = F)

#table <- read.csv2("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112.bed", sep = "\t", header = F)
#table[table$V5 == "","V5"] <- "NA"
#write.table(table, "/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_biomart_protein_coding_genes_112_filled.bed", quote = F, sep = "\t",row.names = F, col.names = F)





data <- read.delim('quantification_deeptools/quantification_edu_genes_plus1kb_fig4f_allgenes.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(data)[1] <- "chr"


data$edu_B10_1 <- data$edu_B10_1 * 1e9 / ((data$end -data$start) * as.numeric(system("grep edu_B10_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$edu_B10_2 <- data$edu_B10_2 * 1e9 / ((data$end -data$start) * as.numeric(system("grep edu_B10_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$edu_B12_1 <- data$edu_B12_1 * 1e9 / ((data$end -data$start) * as.numeric(system("grep edu_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$edu_B12_2 <- data$edu_B12_2 * 1e9 / ((data$end -data$start) * as.numeric(system("grep edu_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$edu_B2_1 <- data$edu_B2_1 * 1e9 / ((data$end -data$start) * as.numeric(system("grep edu_B2_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$edu_B2_2 <- data$edu_B2_2 * 1e9 / ((data$end -data$start) * as.numeric(system("grep edu_B2_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))


data$K37R <- rowMeans(cbind(data$edu_B10_1,data$edu_B10_2))
data$WT <- rowMeans(cbind(data$edu_B12_1,data$edu_B12_2))
data$K37R2  <- rowMeans(cbind(data$edu_B2_1,data$edu_B2_2))

data[,c("edu_B10_1","edu_B10_2","edu_B12_1","edu_B12_2",
        "edu_B2_1","edu_B2_2")] <- NULL

summary(data$WT)

data <- unique(left_join(data, gene.bodies.1kb[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many"))

data <- data[rowSums(data[,4:5]) >0,]



data[3:6] <- log2(data[3:6]+1)

data$quartile <- as.factor(quant_groups(data$WT, ngroups = 4))
data$quartile <- "Low"
data$quartile[data$WT >= 0.5] <- "Medium"
#data$quartile[data$WT >= 1.5] <- "mid"
data$quartile[data$WT >= 2] <- "High"
data$quartile <- factor(data$quartile, levels = c("Low","Medium","High"))
table(data$quartile)


wilcox.test(data$WT[data$quartile == "Low"],data$K37R[data$quartile == "Low"], paired = T)
wilcox.test(data$WT[data$quartile == "Medium"],data$K37R[data$quartile == "Medium"], paired = T)
wilcox.test(data$WT[data$quartile == "High"],data$K37R[data$quartile == "High"], paired = T)


long_data <- data %>% pivot_longer(cols = c(K37R, WT), names_to = "Genotype", values_to = "RPKM")
long_data$Genotype <- factor(long_data$Genotype, levels = c("WT", "K37R"))


# Crear el gráfico
pdf("plots/boxplots/mutants/edu_B12_B10_allgenes_plus1kb_fig4f.pdf", height = 4, width = 6)
ggplot(long_data, aes(x = quartile, y = RPKM, fill = Genotype)) +
  geom_boxplot(outliers = F, staplewidth = 0.5) +  # Quitar los outliers si no los quieres
  stat_compare_means(aes(group = Genotype),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 4.5,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) +  # Ajusta el tamaño del texto del p-valor
  labs(y="EdU-seq-HU signal (RPKM)", x = "") +  
  scale_fill_manual(values = c("#A5272B", "#D19883")) +  # Colores personalizados
  theme_bw() +  # Tema blanco y negro
  theme(text = element_text(size = 12),  # Tamaño de texto general
        axis.text.x = element_text(angle = 45, hjust = 1)) +
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
    
  )  
dev.off()





up <- read.csv2('genes/fc_scatter/groups_allgenes_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_allgenes_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")

data$gene_type <- "ns"
data$gene_type[data$gene_id %in% up$gene_id] <- "Group I"
data$gene_type[data$gene_id %in% down$gene_id] <- "Group II"


df_long <- data %>%
  filter(gene_type %in% c("Group I","Group II")) %>%
  pivot_longer(
    cols = c(WT),
    names_to = "genotype",
    values_to = "FC"
  ) 

df_long$gene_type <- factor(df_long$gene_type, levels = c("WT", "H3.3K37R"))

pdf("plots/boxplots/groups_edu/edu_B12_groups_fc1.pdf", width = 4, height = 7)
ggplot(df_long, aes(x = gene_type, y = FC, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("WT" = "#A5272B",
                               "H3.3K37R" = "#D19883")) +
  stat_compare_means(aes(group = gene_type),
                     method = "wilcox.test",   # O "t.test" si prefieres
                     label = "p.format",
                     label.y = 4,  # Ajusta la posición de los p-value
                     size = 4) +  # Ajusta el tamaño del texto
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "EdU-seq-HU",
    fill = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + coord_cartesian(ylim = c(0, 5))
dev.off()








############################## H3K37me1 BOXPLOTS AND SCATTERPLOTS ##############################






h3k37me1 <- read.delim('quantification_deeptools/quantification_h3k37me1_mut_allgenes.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(h3k37me1)[c(1,4,5)] <- c("chr", "WT", "K37R")


sum(rowSums(h3k37me1[,4:5] > 0) <1)

h3k37me1 <- h3k37me1[rowSums(h3k37me1[,4:5] > 0) >0,]

h3k37me1 <- left_join(h3k37me1, gene.bodies[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")




up <- read.csv2('genes/fc_scatter/groups_allgenes_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_allgenes_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")




h3k37me1$gene_type <- "ns"
h3k37me1$gene_type[h3k37me1$gene_id %in% up$gene_id] <- "up"
h3k37me1$gene_type[h3k37me1$gene_id %in% down$gene_id] <- "down"










df_long <- h3k37me1 %>%
  filter(gene_type %in% c("up", "down")) %>%
  pivot_longer(
    cols = c(WT, K37R),
    names_to = "genotype",
    values_to = "Signal"
  ) %>%
  mutate(
    genotype = recode(genotype,
                      WT_FC = "WT",
                      K37R2_FC = "K37R"),
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


# Crear el gráfico
pdf("plots/boxplots/mutants/h3k37me1_B12_B10_groups_allgenes_fc1.pdf", height = 4, width = 6)
ggplot(df_long, aes(x = gene_type, y = Signal, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("WT" = "#A5272B",
                               "K37R" = "#D19883")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "H3K37me1 signal (RPGC)",
    fill = "Group"
  ) + ylim(c(0, 8.5))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
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
                        label.y = 8,  # Ajusta la posición del p-valor en el gráfico
                        size = 4) 
dev.off()



long_data <- h3k37me1 %>% filter(gene_type %in% "up")


pdf(paste0("plots/scatterplots/scatterplot_h3k37me1_mut_group_down_allgenes_fc1.pdf"))
ggplot(long_data, aes(x = WT, y = K37R)) +
  geom_pointdensity(size = 1, adjust = 4, method = "kde2d") +
  geom_abline(slope = 1, intercept = 0, color ="black", linetype="dashed",linewidth = 1) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data$WT), y = max(data$K37R),
           label = paste0("r = ", round(cor(data$WT, data$K37R, method = "pearson"), 2)),
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


############################## ChRNA BOXPLOTS AND SCATTERPLOTS ##############################


rna <- read.delim('quantification_deeptools/quantification_chrna_B12_B10_allgenes_biomart.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(rna)[1] <- "chr"
rna$width <- rna$end -rna$start
rna$chrna_B12_1 <- rna$chrna_B12_1 * 1e9 / (rna$width * as.numeric(system("grep chrna_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$chrna_B12_2 <- rna$chrna_B12_2 * 1e9 / (rna$width * as.numeric(system("grep chrna_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$chrna_B10_1 <- rna$chrna_B10_1 * 1e9 / (rna$width * as.numeric(system("grep chrna_B10_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$chrna_B10_2 <- rna$chrna_B10_2 * 1e9 / (rna$width * as.numeric(system("grep chrna_B10_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$WT   <- rowMeans(rna[,6:7])
rna$H3.3K37R <- rowMeans(rna[,4:5])
rna <- unique(left_join(rna[,-c(4:7)], gene.bodies[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many"))

# Matrix filtering
rna <- rna[rowSums(rna[,5:6]) > 0,]
rna <- rna[rna$width >= 1000,]



up <- read.csv2('genes/fc_scatter/groups_allgenes_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_allgenes_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")




rna$gene_type <- "ns"
rna$gene_type[rna$gene_id %in% up$gene_id] <- "up"
rna$gene_type[rna$gene_id %in% down$gene_id] <- "down"










df_long <- rna %>%
  filter(gene_type %in% c("up", "down")) %>%
  pivot_longer(
    cols = c(WT, H3.3K37R),
    names_to = "genotype",
    values_to = "Signal"
  ) %>%
  mutate(
    genotype = recode(genotype,
                      WT_FC = "WT",
                      H3.3K37R2_FC = "H3.3K37R"),
    group = case_when(
      genotype == "WT" & gene_type == "up" ~ "WT - UP",
      genotype == "H3.3K37R" & gene_type == "up" ~ "H3.3K37R - UP",
      genotype == "WT" & gene_type == "down" ~ "WT - DOWN",
      genotype == "H3.3K37R" & gene_type == "down" ~ "H3.3K37R - DOWN"
    ),
    group = factor(group, levels = c(
      "WT - UP",
      "H3.3K37R - UP",
      "WT - DOWN",
      "H3.3K37R - DOWN"
    ))
  )


# Crear el gráfico
pdf("plots/boxplots/mutants/rna_B12_B10_groups_allgenes_fc1.pdf", height = 4, width = 6)
ggplot(df_long, aes(x = group, y = Signal, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("WT" = "#A5272B",
                               "H3.3K37R" = "#D19883")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "ChRNA signal (RPKM)",
    fill = "Group"
  ) + ylim(c(0, 4.5))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
    
  )  
  stat_compare_means(aes(group = genotype),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 8.2,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) 
dev.off()



long_data <- rna %>% filter(gene_type %in% "up")


pdf(paste0("plots/scatterplots/scatterplot_rna_mut_group_down_allgenes_fc1.pdf"))
ggplot(long_data, aes(x = WT, y = H3.3K37R)) +
  geom_pointdensity(size = 1, adjust = 4, method = "kde2d") +
  geom_abline(slope = 1, intercept = 0, color ="black", linetype="dashed",linewidth = 1) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data$WT), y = max(data$H3.3K37R),
           label = paste0("r = ", round(cor(data$WT, data$H3.3K37R, method = "pearson"), 2)),
           size = 5, fontface = "bold", color = "#de1f09") +
  theme_minimal(base_size = 14) +
  labs(  x = "WT (log2 RPGC)",
         y = "H3.3H3.3K37R (log2 RPGC)") +
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



df_long <- rna %>%
  pivot_longer(
    cols = c(WT, H3.3K37R),
    names_to = "genotype",
    values_to = "Signal"
  ) %>% mutate(genotype = factor(genotype, levels = c("WT", "H3.3K37R")))


pdf("plots/boxplots/mutants/chrna_B12_B10_allgenes_biomart.pdf", height = 4, width = 4)
ggplot(df_long, aes(x = genotype, y = Signal, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("WT" = "#A5272B",
                               "H3.3K37R" = "#D19883")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "ChRNA-seq",
    fill = "Group"
  ) + coord_cartesian(ylim = c(0, 4.5))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
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
                   label.y = 4,  # Ajusta la posición del p-valor en el gráfico
                   size = 4) 
dev.off()


############################## MCM-2 BOXPLOTS AND SCATTERPLOTS ##############################


mcm2 <- read.delim('quantification_deeptools/quantification_mcm2_mut_allgenes.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(mcm2)[c(1,4,5)] <- c("chr", "WT", "K37R")


sum(rowSums(mcm2[,4:5] > 0) <1)

mcm2 <- mcm2[rowSums(mcm2[,4:5] > 0) >0,]

mcm2 <- left_join(mcm2, gene.bodies[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")



rna <- read.delim('quantification_deeptools/chrna_B12_genebodies_rawcounts.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(rna)[1] <- "chr"
rna$width <- rna$end -rna$start
rna$chrna_B12_1 <- rna$chrna_B12_1 * 1e9 / (rna$width * as.numeric(system("grep chrna_B12_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$chrna_B12_2 <- rna$chrna_B12_2 * 1e9 / (rna$width * as.numeric(system("grep chrna_B12_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
rna$rna <- rowMeans(rna[,4:5])
rna$chr <- as.numeric(rna$chr)



rna <- left_join(rna[,-c(4:6)], gene.bodies[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many")


data <- left_join(mcm2[,4:6],rna[,4:5], by="gene_id", relationship = "many-to-many")
data <- data[rowSums(data[,c(1,2,4)]) > 0,]
#data[,c(1,2,4)] <- log2(data[,c(1,2,4)]+1)


data$quartile <- as.factor(quant_groups(data$rna, ngroups = 4))

summary(data$rna)
data$quartile <- "no"
data$quartile[data$rna >= 0.1] <- "low"
data$quartile[data$rna >= 0.5] <- "mid"
data$quartile[data$rna >= 1.5] <- "high"
data$quartile <- factor(data$quartile, levels = c("no","low","mid","high"))
table(data$quartile)












long_mcm2 <- data %>% pivot_longer(cols = c(K37R, WT), names_to = "Genotype", values_to = "RPKM")%>%
  filter(Genotype %in% c("WT")) %>%
  filter(quartile %in% c("Q1", "Q2", "Q3", "Q4"))
long_mcm2$Genotype <- factor(long_mcm2$Genotype, levels = c("WT", "K37R"))



# Crear el gráfico
pdf("plots/boxplots/mutants/mcm2_B12_quartiles_rna.pdf", height = 4, width = 6)
ggplot(long_mcm2, aes(x = quartile, y = RPKM, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +  # Quitar los outliers si no los quieres
  stat_compare_means(aes(group = quartile),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 6,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) +  # Ajusta el tamaño del texto del p-valor
  labs(y="MCM-2 signal (RPGC)", x = "") +  
  scale_fill_manual(values = c("#A5272B", "#D19883")) +  # Colores personalizados
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + coord_cartesian(ylim = c(0, 5))+
  coord_cartesian(ylim = c(0.2,2.5))
dev.off()


pdf("plots/boxplots/mutants/mcm2_B12_B10_allgenes.pdf", height = 4, width = 3)
ggplot(long_mcm2, aes(x = Genotype, y = RPKM, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +  # Quitar los outliers si no los quieres
  stat_compare_means(aes(group = Genotype),  # Comparar entre los dos genotipos
                     method = "wilcox.test",  # O "t.test" si prefieres
                     label = "p.format",  # Mostrar el p-valor
                     label.y = 2.2,  # Ajusta la posición del p-valor en el gráfico
                     size = 4) +  # Ajusta el tamaño del texto del p-valor
  labs(y="MCM-2 signal (RPGC)", x = "") +  
  scale_fill_manual(values = c("#A5272B", "#D19883")) +  # Colores personalizados
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    
    # 🔹 Quitar grid vertical
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + coord_cartesian(ylim = c(0, 5))+
  coord_cartesian(ylim = c(0.2,2.5))
dev.off()









up <- read.csv2('genes/fc_scatter/groups_allgenes_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_allgenes_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")










mcm2$gene_type <- "ns"
mcm2$gene_type[mcm2$gene_id %in% up$gene_id] <- "Group I"
mcm2$gene_type[mcm2$gene_id %in% down$gene_id] <- "Group II"










df_long <- mcm2 %>%
  filter(gene_type %in% c("Group I", "Group II")) %>%
  pivot_longer(
    cols = c(WT, K37R),
    names_to = "genotype",
    values_to = "Signal"
  ) 
df_long$genotype <- factor(df_long$genotype, levels = c("WT", "K37R"))



# Crear el gráfico
pdf("plots/boxplots/mutants/mcm2_B12_B10_groups_allgenes_fc1.pdf", height = 4, width = 5)
ggplot(df_long, aes(x = gene_type, y = Signal, fill = genotype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("WT" = "#A5272B",
                               "K37R" = "#D19883")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "MCM-2 signal (RPGC)",
    fill = "Group"
  ) + coord_cartesian(ylim=c(0, 2.2))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 🔹 Líneas horizontales suaves
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.5),
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



long_data <- mcm2 %>% filter(gene_type %in% "up")


pdf(paste0("plots/scatterplots/scatterplot_mcm2_mut_group_down_allgenes_fc1.pdf"))
ggplot(long_data, aes(x = WT, y = K37R)) +
  geom_pointdensity(size = 1, adjust = 4, method = "kde2d") +
  geom_abline(slope = 1, intercept = 0, color ="black", linetype="dashed",linewidth = 1) +
  geom_smooth(method = "lm", color = "#de1f09", se = FALSE, linewidth = 0.8) +
  annotate("text",
           x = max(data$WT), y = max(data$K37R),
           label = paste0("r = ", round(cor(data$WT, data$K37R, method = "pearson"), 2)),
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




############################## REPLI-SEQ S1-4 BOXPLOTS ##############################




data <- read.delim('quantification_deeptools/quantification_repli_allgenes.txt', sep = "\t", quote = "'",  header = T,stringsAsFactors = T )
colnames(data)[1] <- "chr"
data$width <- data$end -data$start
data$repli_s1_1 <- data$repli_s1_1 * 1e9 / (data$width * as.numeric(system("grep repli_s1_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$repli_s1_2 <- data$repli_s1_2 * 1e9 / (data$width * as.numeric(system("grep repli_s1_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$repli_s2_1 <- data$repli_s2_1 * 1e9 / (data$width * as.numeric(system("grep repli_s2_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$repli_s2_2 <- data$repli_s2_2 * 1e9 / (data$width * as.numeric(system("grep repli_s2_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$repli_s3_1 <- data$repli_s3_1 * 1e9 / (data$width * as.numeric(system("grep repli_s3_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$repli_s3_2 <- data$repli_s3_2 * 1e9 / (data$width * as.numeric(system("grep repli_s3_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$repli_s4_1 <- data$repli_s4_1 * 1e9 / (data$width * as.numeric(system("grep repli_s4_1 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$repli_s4_2 <- data$repli_s4_2 * 1e9 / (data$width * as.numeric(system("grep repli_s4_2 ../bam_files/per_replicate/total_counts.tsv | cut -f2",intern = T)))
data$S1 <- rowMeans(data[,4:5])
data$S2 <- rowMeans(data[,6:7])
data$S3 <- rowMeans(data[,8:9])
data$S4 <- rowMeans(data[,10:11])

data <- unique(left_join(data[,-c(4:11)], gene.bodies[,1:4],  by = c("chr", "start", "end"), relationship = "many-to-many"))

# Matrix filtering
data <- data[rowSums(data[,5:8]) > 0,]
#data <- data[data$width >= 1000,]



up <- read.csv2('genes/fc_scatter/groups_allgenes_B10_up_fc1.bed', sep = '\t', header = F)
down <- read.csv2('genes/fc_scatter/groups_allgenes_B10_down_fc1.bed', sep = '\t', header = F)
colnames(up) <- colnames(down) <- c("chr", "start","end","gene_id","score","strand")




data$Group <- "ns"
data$Group[data$gene_id %in% up$gene_id] <- "Group I"
data$Group[data$gene_id %in% down$gene_id] <- "Group II"










df_long <- data %>%
  filter(Group %in% c("Group I", "Group II")) %>%
  pivot_longer(
    cols = c(S1, S2, S3, S4),
    names_to = "Phase",
    values_to = "Signal"
  ) 


# Crear el gráfico
pdf("plots/boxplots/groups_edu/repli_groups_allgenes_fc1.pdf", height = 4, width = 6)
ggplot(df_long, aes(x = Phase, y = Signal, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("Group I" = "#E74C3C", "Group II"= "#3498DB")) +
  #scale_y_continuous(    limits = c(-3, 1.25),    breaks = seq(-3, 1.25, by = 0.75)  ) +
  labs(
    x = "",
    y = "Repli-seq",
    fill = "Group"
  ) + coord_cartesian(ylim=c(0, 2.8))+
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # 🔹 Línea del eje Y
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    
    # 🔹 Activar ticks del eje Y
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm")  ) +
stat_compare_means(aes(group = Group),  # Comparar entre los dos genotipos
                   method = "wilcox.test",  # O "t.test" si prefieres
                   label = "p.format",  # Mostrar el p-valor
                   label.y = 2.6,  # Ajusta la posición del p-valor en el gráfico
                   size = 4) 
dev.off()




