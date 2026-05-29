
############################## LOADING DATA AND FUNCTIONS ##############################

# Setting working directory and loading all necessary packages and data.

setwd('/media/gmzlab/4TB/projects/h3k37/results/')
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
#suppressMessages(library(gUtils))
#suppressMessages(library(tidyverse))
#BiocManager::install("ChIPseeker")
suppressMessages(library(ChIPseeker))
#suppressMessages(library(ggupset))
#suppressMessages(library(GenomicRanges))
#suppressMessages(library(ggpubr))
#suppressMessages(library(cowplot))
#suppressMessages(library(heatmaply))
suppressMessages(library(ggplot2))
#BiocManager::install("txdbmaker")
#suppressMessages(library(txdbmaker))
#install.packages("eulerr")
suppressMessages(library(eulerr))
suppressMessages(library(bedtoolsr))
suppressMessages(library(org.Hs.eg.db))
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
#suppressMessages(library(GenomicFeatures))
#BiocManager::install("EnsDb.Hsapiens.v86")
suppressMessages(library(EnsDb.Hsapiens.v86))


gene.id <- mapIds(
  org.Hs.eg.db,
  keys = read.csv2("genes/genes_chromHMM/genes_broad_E5.bed", sep = "\t", header = F)[,4],
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
sum(is.na(gene.id))# eliminar NA (genes sin conversión)
gene.id <- gene.id[!is.na(gene.id)]



up.gene.id <- mapIds(
  org.Hs.eg.db,
  keys = read.csv2("genes/fc_scatter/groups_B10_up_fc1.bed", sep = "\t", header = F)[,4],
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
sum(is.na(up.gene.id))# eliminar NA (genes sin conversión)
up.gene.id <- up.gene.id[!is.na(up.gene.id)]

down.gene.id <- mapIds(
  org.Hs.eg.db,
  keys = read.csv2("genes/fc_scatter/groups_B10_down_fc1.bed", sep = "\t", header = F)[,4],
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
sum(is.na(down.gene.id))# eliminar NA (genes sin conversión)
down.gene.id <- down.gene.id[!is.na(down.gene.id)]



############################## RARE CHROMOSOMES REMOVAL ##############################

peaks <- readPeakFile("peaks/edu_B12/edu_B12_peak_fc7.narrowPeak")
peaks <- peaks[seqnames(peaks) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")]
seqlevels(peaks) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
rtracklayer::export(peaks, "peaks/chrom_curated/edu_B12_peaks.bed")

peaks <- readPeakFile("peaks/edu_B10/edu_B10_peak_fc7.narrowPeak")
peaks <- peaks[seqnames(peaks) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")]
seqlevels(peaks) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
rtracklayer::export(peaks, "peaks/chrom_curated/edu_B10_peaks.bed")

peaks <- readPeakFile("peaks/edu_B2/edu_B2_peak_fc7.narrowPeak")
peaks <- peaks[seqnames(peaks) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")]
seqlevels(peaks) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
rtracklayer::export(peaks, "peaks/chrom_curated/edu_B2_peaks.bed")

peaks <- readPeakFile("/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed")
peaks <- peaks[seqnames(peaks) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")]
seqlevels(peaks) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
rtracklayer::export(peaks, "peaks/chrom_curated/Homo_sapiens_ann_genes_112.bed")

peaks <- readPeakFile("genes/fc_scatter/groups_B10_up_fc1_plus10kb.bed")
peaks <- peaks[seqnames(peaks) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")]
seqlevels(peaks) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
rtracklayer::export(peaks, "peaks/chrom_curated/groups_B10_up_fc1_plus10kb.bed")

peaks <- readPeakFile("genes/fc_scatter/groups_B10_down_fc1_plus10kb.bed")
peaks <- peaks[seqnames(peaks) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")]
seqlevels(peaks) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
rtracklayer::export(peaks, "peaks/chrom_curated/groups_B10_down_fc1_plus10kb.bed")


############################## INTERSECTIONS ##############################

system("bedtools intersect -u -a peaks/chrom_curated/edu_B12_peaks.bed -b peaks/chrom_curated/edu_B10_peaks.bed > peaks/chrom_curated/edu_B12_B10_peaks.bed  ")
system("bedtools intersect -v -a peaks/chrom_curated/edu_B12_peaks.bed -b peaks/chrom_curated/edu_B10_peaks.bed > peaks/chrom_curated/edu_B12_noB10_peaks.bed")
system("bedtools intersect -u -a peaks/chrom_curated/edu_B10_peaks.bed -b peaks/chrom_curated/edu_B12_peaks.bed > peaks/chrom_curated/edu_B10_B12_peaks.bed  ")
system("bedtools intersect -v -a peaks/chrom_curated/edu_B10_peaks.bed -b peaks/chrom_curated/edu_B12_peaks.bed > peaks/chrom_curated/edu_B10_noB12_peaks.bed")
system("bedtools intersect -u -a peaks/chrom_curated/edu_B12_peaks.bed -b peaks/chrom_curated/edu_B2_peaks.bed  > peaks/chrom_curated/edu_B12_B2_peaks.bed   ")
system("bedtools intersect -v -a peaks/chrom_curated/edu_B12_peaks.bed -b peaks/chrom_curated/edu_B2_peaks.bed  > peaks/chrom_curated/edu_B12_noB2_peaks.bed ")
system("bedtools intersect -u -a peaks/chrom_curated/edu_B2_peaks.bed  -b peaks/chrom_curated/edu_B12_peaks.bed > peaks/chrom_curated/edu_B2_B12_peaks.bed   ")
system("bedtools intersect -v -a peaks/chrom_curated/edu_B2_peaks.bed  -b peaks/chrom_curated/edu_B12_peaks.bed > peaks/chrom_curated/edu_B2_noB12_peaks.bed ")



############################## PEAKS LOADING ##############################

B12.peaks       <- readPeakFile('peaks/chrom_curated/edu_B12_peaks.bed'      )
B12.B10.peaks   <- readPeakFile('peaks/chrom_curated/edu_B12_B10_peaks.bed'  )
B12.B2.peaks    <- readPeakFile('peaks/chrom_curated/edu_B12_B2_peaks.bed'   )
B12.noB10.peaks <- readPeakFile('peaks/chrom_curated/edu_B12_noB10_peaks.bed')
B12.noB2.peaks  <- readPeakFile('peaks/chrom_curated/edu_B12_noB2_peaks.bed' )

B10.peaks       <- readPeakFile('peaks/chrom_curated/edu_B10_peaks.bed'      )
B10.B12.peaks   <- readPeakFile('peaks/chrom_curated/edu_B10_B12_peaks.bed'  )
B10.noB12.peaks <- readPeakFile('peaks/chrom_curated/edu_B10_noB12_peaks.bed')

B2.peaks        <- readPeakFile('peaks/chrom_curated/edu_B2_peaks.bed'       )
B2.B12.peaks    <- readPeakFile('peaks/chrom_curated/edu_B2_B12_peaks.bed'   )
B2.noB12.peaks  <- readPeakFile('peaks/chrom_curated/edu_B2_noB12_peaks.bed' )

############################## PEAKS SELECTION ACCORDING TO WIDTH ##############################


B12.peaks$width       <- width(B12.peaks      )
B12.B10.peaks$width   <- width(B12.B10.peaks  )
B12.B2.peaks$width    <- width(B12.B2.peaks   )
B12.noB10.peaks$width <- width(B12.noB10.peaks)
B12.noB2.peaks$width  <- width(B12.noB2.peaks )
B10.peaks$width       <- width(B10.peaks      )
B10.B12.peaks$width   <- width(B10.B12.peaks  )
B10.noB12.peaks$width <- width(B10.noB12.peaks)
B2.peaks$width        <- width(B2.peaks       )
B2.B12.peaks$width    <- width(B2.B12.peaks   )
B2.noB12.peaks$width  <- width(B2.noB12.peaks )



length(B12.peaks[B12.peaks$width >= 2000])
length(B12.peaks[B12.peaks$width <  2000])

length(B12.noB10.peaks[B12.noB10.peaks$width >= 2000])
length(B12.noB10.peaks[B12.noB10.peaks$width <  2000])



B12.peaks <- readPeakFile('peaks/chrom_curated/edu_B12_peaks.bed')
B12.peaks <- readPeakFile('peaks/chrom_curated/edu_B12_peaks.bed')
B10.peaks <- readPeakFile('peaks/chrom_curated/edu_B10_peaks.bed')
B10.peaks <- readPeakFile('peaks/chrom_curated/edu_B10_peaks.bed')
B2.peaks  <- readPeakFile('peaks/chrom_curated/edu_B2_peaks.bed' )
B2.peaks  <- readPeakFile('peaks/chrom_curated/edu_B2_peaks.bed' )



############################## PEAKS ANNOTATION AND PIECHART (ONE PEAKS FILE) ##############################

#peak.file <- readPeakFile('peaks/chrom_curated/edu_B2_noB12.bed')

peak.file <- B12.peaks
pdf.file  <- 'plots/barplots/edu_B12_peaks_annotation.pdf'

seqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene)
#seqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene) <- sub("^chr", "", seqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene))
#seqlevels(readPeakFile('peaks/chrom_curated/edu_B2_noB12_peaks.bed'))
df <- as.data.frame(annotatePeak( B12.peaks, tssRegion=c(-1000, 1000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene))
df$region <- sub(" \\(.*", "", df$annotation)
df$group <- ifelse(df$region %in%  c("Intron", "Exon", "3' UTR", "5' UTR"), "Gene body",
            ifelse(df$region %in% c("Distal Intergenic", "Downstream"), "Intergenic", df$region))
df_plot <- df %>% dplyr::count(group) %>% mutate(porcentaje = round(100 * n / sum(n), 1),
           etiqueta = paste0(group, "\n", porcentaje, "%"))
pdf(pdf.file)
ggplot(df_plot, aes(x = "", y = porcentaje, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  geom_text(aes(label = etiqueta), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("#429A48",  "#e39610",  "#4f95ff" )) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

#plotDistToTSS(annotatePeak(peak.file, tssRegion=c(-1000, 1000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene))


############################## PEAKS ANNOTATION AND BARPLOT (TWO PEAKS FILES) ##############################

process_peaks <- function(peaks) {
  df <- as.data.frame(
    annotatePeak(peaks, 
                 tssRegion = c(-1000, 0), 
                 TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
  )
  
  # igual que en tu código
  df$region <- sub(" \\(.*", "", df$annotation)
  
  df$group <- ifelse(
    df$region %in% c("Intron", "Exon", "3' UTR", "5' UTR"), "Genic",
    ifelse(df$region %in% c("Distal Intergenic", "Downstream", "Promoter"), 
           "Intergenic",
           df$region)
  )
  
  df %>% count(group) %>% mutate(percent = 100 * n / sum(n))
}

# ------------ B10
wt_df  <- mutate(process_peaks(B12.B10.peaks), condition ="WT")
mut_df <- mutate(process_peaks(B10.noB12.peaks), condition ="H3.3K37R")
df_all <- rbind(wt_df, mut_df)
#df_all$group <- factor(df_all$group, levels = c("Gene body","Promoter","Intergenic"))
df_all$group <- factor(df_all$group, levels = c("Genic","Intergenic"))
df_all$condition <- factor(df_all$condition, levels=c("WT","H3.3K37R"))

pdf("plots/barplots/annotation_B10noB12_vs_B12andB10_peaks.pdf",height = 8, width = 4)
ggplot(df_all, aes(x=condition, y=percent, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Genic"="#e39610", "Intergenic"="#aaaaaa")) +
  xlab("") +
  ylab("% peaks") +
  ylim(0,77)+
  theme_bw(base_size=14) +
  #ggtitle("Distribution of peak annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank())
dev.off()

write.table(df_all, "edu_peaks_annotation_B10_B12.txt", quote = F, sep = "\t", row.names = F)

# ------------ B2
wt_df  <- mutate(process_peaks(B12.B2.peaks), condition ="WT")
mut_df <- mutate(process_peaks(B2.noB12.peaks), condition ="H3.3K37R #2")
df_all <- rbind(wt_df, mut_df)
#df_all$group <- factor(df_all$group, levels = c("Gene body","Promoter","Intergenic"))
df_all$group <- factor(df_all$group, levels = c("Genic","Intergenic"))
df_all$condition <- factor(df_all$condition, levels=c("WT","H3.3K37R #2"))

pdf("plots/barplots/annotation_B2noB12_vs_B12andB2_peaks.pdf",height = 8, width = 4)
ggplot(df_all, aes(x=condition, y=percent, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Genic"="#e39610", "Intergenic"="#aaaaaa")) +
  xlab("") +
  ylab("% peaks") +
  ylim(0,77)+
  theme_bw(base_size=14) +
  #ggtitle("Distribution of peak annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank())
dev.off()


write.table(df_all, "edu_peaks_annotation_B2_B12.txt", quote = F, sep = "\t", row.names = F)



############################## PEAKS ANNOTATION IN EUCHROMATIN AND BARPLOT (TWO PEAKS FILES) ##############################

process_peaks <- function(peaks) {
  
  # --- Anotación normal con TxDb completo ---
  df <- as.data.frame(
    annotatePeak(
      peaks, 
      tssRegion = c(-1000, 0), 
      TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
    )
  )
  
  # --- FILTRO SOLO A TUS GENES EUCROMÁTICOS (OPCIÓN A) ---
  df <- df[df$geneId %in% gene.id, ]
  
  # --- Tu procesamiento original ---
  df$region <- sub(" \\(.*", "", df$annotation)
  
  df$group <- ifelse(
    df$region %in% c("Intron", "Exon", "3' UTR", "5' UTR"), "Genic",
    ifelse(df$region %in% c("Distal Intergenic", "Downstream", "Promoter"), 
           "Intergenic",
           df$region)
  )
  
  df %>% count(group) %>% mutate(percent = 100 * n / sum(n))
}


# ------------ B10
wt_df  <- mutate(process_peaks(B12.B10.peaks), condition ="WT")
mut_df <- mutate(process_peaks(B10.noB12.peaks), condition ="H3.3K37R")
df_all <- rbind(wt_df, mut_df)
#df_all$group <- factor(df_all$group, levels = c("Gene body","Promoter","Intergenic"))
df_all$group <- factor(df_all$group, levels = c("Genic","Intergenic"))
df_all$condition <- factor(df_all$condition, levels=c("WT","H3.3K37R"))

pdf("plots/barplots/annotation_B10noB12_vs_B12andB10_peaks.pdf",height = 8, width = 4)
ggplot(df_all, aes(x=condition, y=percent, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Genic"="#e39610", "Intergenic"="#aaaaaa")) +
  xlab("") +
  ylab("% peaks") +
  ylim(0,85)+
  theme_bw(base_size=14) +
  #ggtitle("Distribution of peak annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank())
dev.off()

write.table(df_all, "edu_peaks_annotation_B10_B12.txt", quote = F, sep = "\t", row.names = F)

# ------------ B2
wt_df  <- mutate(process_peaks(B12.B2.peaks), condition ="WT")
mut_df <- mutate(process_peaks(B2.noB12.peaks), condition ="H3.3K37R #2")
df_all <- rbind(wt_df, mut_df)
#df_all$group <- factor(df_all$group, levels = c("Gene body","Promoter","Intergenic"))
df_all$group <- factor(df_all$group, levels = c("Genic","Intergenic"))
df_all$condition <- factor(df_all$condition, levels=c("WT","H3.3K37R #2"))

pdf("plots/barplots/annotation_B2noB12_vs_B12andB2_peaks.pdf",height = 8, width = 4)
ggplot(df_all, aes(x=condition, y=percent, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Genic"="#e39610", "Intergenic"="#aaaaaa")) +
  xlab("") +
  ylab("% peaks") +
  ylim(0,85)+
  theme_bw(base_size=14) +
  #ggtitle("Distribution of peak annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank())
dev.off()


write.table(df_all, "edu_peaks_annotation_B2_B12.txt", quote = F, sep = "\t", row.names = F)









############################## PEAKS ANNOTATION IN EDU GROUP UP AND BARPLOT (TWO PEAKS FILES) ##############################

process_peaks <- function(peaks) {
  
  # --- Anotación normal con TxDb completo ---
  df <- as.data.frame(
    annotatePeak(
      peaks, 
      tssRegion = c(-1000, 0), 
      TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
    )
  )
  
  # --- FILTRO SOLO A TUS GENES EUCROMÁTICOS (OPCIÓN A) ---
  df <- df[df$geneId %in% up.gene.id, ]
  
  # --- Tu procesamiento original ---
  df$region <- sub(" \\(.*", "", df$annotation)
  
  df$group <- ifelse(
    df$region %in% c("Intron", "Exon", "3' UTR", "5' UTR"), "Genic",
    ifelse(df$region %in% c("Distal Intergenic", "Downstream", "Promoter"), 
           "Intergenic",
           df$region)
  )
  
  df %>% count(group) %>% mutate(percent = 100 * n / sum(n))
}


# ------------ B10
wt_df  <- mutate(process_peaks(B12.B10.peaks), condition ="WT")
mut_df <- mutate(process_peaks(B10.noB12.peaks), condition ="H3.3K37R")
df_all <- rbind(wt_df, mut_df)
#df_all$group <- factor(df_all$group, levels = c("Gene body","Promoter","Intergenic"))
df_all$group <- factor(df_all$group, levels = c("Genic","Intergenic"))
df_all$condition <- factor(df_all$condition, levels=c("WT","H3.3K37R"))

pdf("plots/barplots/annotation_B10noB12_vs_B12andB10_peaks_edu_groups_B10_up_fc1.pdf",height = 8, width = 4)
ggplot(df_all, aes(x=condition, y=percent, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Genic"="#e39610", "Intergenic"="#aaaaaa")) +
  xlab("") +
  ylab("% peaks") +
  ylim(0,85)+
  theme_bw(base_size=14) +
  #ggtitle("Distribution of peak annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank())
dev.off()

write.table(df_all, "edu_peaks_annotation_B10noB12_vs_B12andB10_groups_B10_up_fc1.txt", quote = F, sep = "\t", row.names = F)

# ------------ B2
wt_df  <- mutate(process_peaks(B12.B2.peaks), condition ="WT")
mut_df <- mutate(process_peaks(B2.noB12.peaks), condition ="H3.3K37R #2")
df_all <- rbind(wt_df, mut_df)
#df_all$group <- factor(df_all$group, levels = c("Gene body","Promoter","Intergenic"))
df_all$group <- factor(df_all$group, levels = c("Genic","Intergenic"))
df_all$condition <- factor(df_all$condition, levels=c("WT","H3.3K37R #2"))

pdf("plots/barplots/annotation_B2noB12_vs_B12andB2_peaks_groups_B10_up_fc1.pdf",height = 8, width = 4)
ggplot(df_all, aes(x=condition, y=percent, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Genic"="#e39610", "Intergenic"="#aaaaaa")) +
  xlab("") +
  ylab("% peaks") +
  ylim(0,85)+
  theme_bw(base_size=14) +
  #ggtitle("Distribution of peak annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank())
dev.off()


write.table(df_all, "edu_peaks_annotation_B2noB12_vs_B12andB2_edu_groups_B10_up_fc1.txt", quote = F, sep = "\t", row.names = F)









############################## PEAKS ANNOTATION IN EDU GROUP DOWN AND BARPLOT (TWO PEAKS FILES) ##############################

process_peaks <- function(peaks) {
  
  # --- Anotación normal con TxDb completo ---
  df <- as.data.frame(
    annotatePeak(
      peaks, 
      tssRegion = c(-1000, 0), 
      TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene
    )
  )
  
  # --- FILTRO SOLO A TUS GENES EUCROMÁTICOS (OPCIÓN A) ---
  df <- df[df$geneId %in% down.gene.id, ]
  
  # --- Tu procesamiento original ---
  df$region <- sub(" \\(.*", "", df$annotation)
  
  df$group <- ifelse(
    df$region %in% c("Intron", "Exon", "3' UTR", "5' UTR"), "Genic",
    ifelse(df$region %in% c("Distal Intergenic", "Downstream", "Promoter"), 
           "Intergenic",
           df$region)
  )
  
  df %>% count(group) %>% mutate(percent = 100 * n / sum(n))
}


# ------------ B10
wt_df  <- mutate(process_peaks(B12.B10.peaks), condition ="WT")
mut_df <- mutate(process_peaks(B10.noB12.peaks), condition ="H3.3K37R")
df_all <- rbind(wt_df, mut_df)
#df_all$group <- factor(df_all$group, levels = c("Gene body","Promoter","Intergenic"))
df_all$group <- factor(df_all$group, levels = c("Genic","Intergenic"))
df_all$condition <- factor(df_all$condition, levels=c("WT","H3.3K37R"))

pdf("plots/barplots/annotation_B10noB12_vs_B12andB10_peaks_edu_groups_B10_down_fc1.pdf",height = 8, width = 4)
ggplot(df_all, aes(x=condition, y=percent, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Genic"="#e39610", "Intergenic"="#aaaaaa")) +
  xlab("") +
  ylab("% peaks") +
  ylim(0,85)+
  theme_bw(base_size=14) +
  #ggtitle("Distribution of peak annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank())
dev.off()

write.table(df_all, "edu_peaks_annotation_B10noB12_vs_B12andB10_groups_B10_down_fc1.txt", quote = F, sep = "\t", row.names = F)

# ------------ B2
wt_df  <- mutate(process_peaks(B12.B2.peaks), condition ="WT")
mut_df <- mutate(process_peaks(B2.noB12.peaks), condition ="H3.3K37R #2")
df_all <- rbind(wt_df, mut_df)
#df_all$group <- factor(df_all$group, levels = c("Gene body","Promoter","Intergenic"))
df_all$group <- factor(df_all$group, levels = c("Genic","Intergenic"))
df_all$condition <- factor(df_all$condition, levels=c("WT","H3.3K37R #2"))

pdf("plots/barplots/annotation_B2noB12_vs_B12andB2_peaks_groups_B10_down_fc1.pdf",height = 8, width = 4)
ggplot(df_all, aes(x=condition, y=percent, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Genic"="#e39610", "Intergenic"="#aaaaaa")) +
  xlab("") +
  ylab("% peaks") +
  ylim(0,85)+
  theme_bw(base_size=14) +
  #ggtitle("Distribution of peak annotations") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank())
dev.off()


write.table(df_all, "edu_peaks_annotation_B2noB12_vs_B12andB2_edu_groups_B10_down_fc1.txt", quote = F, sep = "\t", row.names = F)









############################## PEAKS OVERLAP ##############################


# B10 vs B12

pdf('plots/vendiagrams/venndiagram_edu_B10_B12.pdf')
plot(euler(c("WT" = nrow(bt.intersect(a = B12.peaks, b = B10.peaks, v = T)),
             "H3.3K37R" = nrow(bt.intersect(a = B10.peaks, b = B12.peaks, v = T)),
             "WT&H3.3K37R" = nrow(bt.intersect(a = B12.peaks, b = B10.peaks, u = T)))),
     fills = list(fill = c("#a5272b", "#D19883"), alpha = 0.8),  # colores y transparencia
     edges = list(col = "black", lwd = 1.5),                     # bordes de círculos
     labels = list(font = 1, cex = 1.5), quantities = TRUE)
dev.off()

nrow(bt.intersect(a = B12.peaks, b = B10.peaks, v = T)) # B12 ONLY
nrow(bt.intersect(a = B12.peaks, b = B10.peaks, u = T)) # B12 with B10 overlap
nrow(bt.intersect(a = B10.peaks, b = B12.peaks, u = T)) # B10 with B12 overlap
nrow(bt.intersect(a = B10.peaks, b = B12.peaks, v = T)) # B10 ONLY

pdf('plots/vendiagrams/venndiagram_edu_B10_B12_adjusted.pdf')
plot(euler(c("WT" = nrow(bt.intersect(a = B12.peaks, b = B10.peaks, v = T)),
             "H3.3K37R" = nrow(bt.intersect(a = B10.peaks, b = B12.peaks, v = T))*nrow(bt.intersect(a = B12.peaks, b = B10.peaks, u = T)) / nrow(bt.intersect(a = B10.peaks, b = B12.peaks, u = T)),
             "WT&H3.3K37R" = nrow(bt.intersect(a = B12.peaks, b = B10.peaks, u = T)))),
     fills = list(fill = c("#a5272b", "#D19883"), alpha = 0.8),  # colores y transparencia
     edges = list(col = "black", lwd = 1.5),                     # bordes de círculos
     labels = list(font = 1, cex = 1.5), quantities = TRUE)
dev.off()



# B2 vs B12

pdf('plots/vendiagrams/venndiagram_edu_B2_B12.pdf')
plot(euler(c("WT" = nrow(bt.intersect(a = B12.peaks, b = B2.peaks, v = T)),
             "H3.3K37R #2" = nrow(bt.intersect(a = B2.peaks, b = B12.peaks, v = T)),
             "WT&H3.3K37R #2" = nrow(bt.intersect(a = B12.peaks, b = B2.peaks, u = T)))),
     fills = list(fill = c("#a5272b", "#6da1c2"), alpha = 0.8),  # colores y transparencia
     edges = list(col = "black", lwd = 1.5),                     # bordes de círculos
     labels = list(font = 1, cex = 1.5), quantities = T )
dev.off()

nrow(bt.intersect(a = B12.peaks, b = B2.peaks, v = T)) # B12 ONLY
nrow(bt.intersect(a = B12.peaks, b = B2.peaks, u = T)) # B12 with B2 overlap
nrow(bt.intersect(a = B2.peaks, b = B12.peaks, u = T)) # B2 with B12 overlap
nrow(bt.intersect(a = B2.peaks, b = B12.peaks, v = T)) # B2 ONLY

pdf('plots/vendiagrams/venndiagram_edu_B2_B12_adjusted.pdf')
plot(euler(c("WT" = nrow(bt.intersect(a = B12.peaks, b = B2.peaks, v = T)),
             "H3.3K37R #2" = nrow(bt.intersect(a = B2.peaks, b = B12.peaks, v = T))*nrow(bt.intersect(a = B12.peaks, b = B2.peaks, u = T)) / nrow(bt.intersect(a = B2.peaks, b = B12.peaks, u = T)),
             "WT&H3.3K37R #2" = nrow(bt.intersect(a = B12.peaks, b = B2.peaks, u = T)))),
     fills = list(fill = c("#a5272b", "#6da1c2"), alpha = 0.8),  # colores y transparencia
     edges = list(col = "black", lwd = 1.5),                     # bordes de círculos
     labels = list(font = 1, cex = 1.5), quantities = TRUE)
dev.off()

# B10 vs B2

pdf('plots/vendiagrams/venndiagram_edu_B10_B2.pdf')
plot(euler(c("H3.3K37R" = nrow(bt.intersect(a = B10.peaks, b = B2.peaks, v = T)),
             "H3.3K37R #2" = nrow(bt.intersect(a = B2.peaks, b = B10.peaks, v = T)),
             "H3.3K37R&H3.3K37R #2" = nrow(bt.intersect(a = B10.peaks, b = B2.peaks, u = T)))),
     fills = list(fill = c("#D19883", "#6da1c2"), alpha = 0.8),  # colores y transparencia
     edges = list(col = "black", lwd = 1.5),                     # bordes de círculos
     labels = list(font = 1, cex = 1.5), quantities = TRUE)
dev.off()

nrow(bt.intersect(a = B10.peaks, b = B2.peaks, v = T)) # B10 ONLY
nrow(bt.intersect(a = B10.peaks, b = B2.peaks, u = T)) # B10 with B2 overlap
nrow(bt.intersect(a = B2.peaks, b = B10.peaks, u = T)) # B2 with B10 overlap
nrow(bt.intersect(a = B2.peaks, b = B10.peaks, v = T)) # B2 ONLY

pdf('plots/vendiagrams/venndiagram_edu_B10_B2_adjusted.pdf')
plot(euler(c("H3.3K37R" = nrow(bt.intersect(a = B10.peaks, b = B2.peaks, v = T)),
             "H3.3K37R #2" = nrow(bt.intersect(a = B2.peaks, b = B10.peaks, v = T))*nrow(bt.intersect(a = B10.peaks, b = B2.peaks, u = T)) / nrow(bt.intersect(a = B2.peaks, b = B10.peaks, u = T)),
             "H3.3K37R&H3.3K37R #2" = nrow(bt.intersect(a = B10.peaks, b = B2.peaks, u = T)))),
     fills = list(fill = c("#D19883", "#6da1c2"), alpha = 0.8),  # colores y transparencia
     edges = list(col = "black", lwd = 1.5),                     # bordes de círculos
     labels = list(font = 1, cex = 1.5), quantities = TRUE)
dev.off()




















# B10 NEW vs B2 NEW


pdf('plots/vendiagrams/venndiagram_edu_B10noB12_B2noB12.pdf')
plot(euler(c("H3.3K37R" = nrow(bt.intersect(a = B10.noB12.peaks, b = B2.noB12.peaks, v = T)),
             "H3.3K37R #2" = nrow(bt.intersect(a = B2.noB12.peaks, b = B10.noB12.peaks, v = T)),
             "H3.3K37R&H3.3K37R #2" = nrow(bt.intersect(a = B10.noB12.peaks, b = B2.noB12.peaks, u = T)))),
     fills = list(fill = c("#D19883", "#6da1c2"), alpha = 0.8),  # colores y transparencia
     edges = list(col = "black", lwd = 1.5),                     # bordes de círculos
     labels = list(font = 1, cex = 1.5), quantities = TRUE)
dev.off()

nrow(bt.intersect(a = B10.noB12.peaks, b = B2.noB12.peaks, v = T)) # B10 ONLY
nrow(bt.intersect(a = B10.noB12.peaks, b = B2.noB12.peaks, u = T)) # B10 with B2 overlap
nrow(bt.intersect(a = B2.noB12.peaks, b = B10.noB12.peaks, u = T)) # B2 with B10 overlap
nrow(bt.intersect(a = B2.noB12.peaks, b = B10.noB12.peaks, v = T)) # B2 ONLY

pdf('plots/vendiagrams/venndiagram_edu_B10_noB12_B2_noB12_adjusted.pdf')
plot(euler(c("H3.3K37R" = nrow(bt.intersect(a = B10.noB12.peaks, b = B2.noB12.peaks, v = T)),
             "H3.3K37R #2" = nrow(bt.intersect(a = B2.noB12.peaks, b = B10.noB12.peaks, v = T))*nrow(bt.intersect(a = B10.noB12.peaks, b = B2.noB12.peaks, u = T)) / nrow(bt.intersect(a = B2.noB12.peaks, b = B10.noB12.peaks, u = T)),
             "H3.3K37R&H3.3K37R #2" = nrow(bt.intersect(a = B10.noB12.peaks, b = B2.noB12.peaks, u = T)))),
     fills = list(fill = c("#D19883", "#6da1c2"), alpha = 0.8),  # colores y transparencia
     edges = list(col = "black", lwd = 1.5),                     # bordes de círculos
     labels = list(font = 1, cex = 1.5), quantities = TRUE)
dev.off()



############################## GENES AND PEAKS ##############################


# Count EdU peaks in annotated genes:
system("bedtools intersect -a /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed -b peaks/chrom_curated/edu_B12_peaks.bed -c > peaks/chrom_curated/count_edu_B12_in_genes.bed")
system("bedtools intersect -a /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed -b peaks/chrom_curated/edu_B10_peaks.bed -c > peaks/chrom_curated/count_edu_B10_in_genes.bed")
system("bedtools intersect -a /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed -b peaks/chrom_curated/edu_B2_peaks.bed  -c > peaks/chrom_curated/count_edu_B2_in_genes.bed")

system("bedtools intersect -a peaks/chrom_curated/groups_B10_up_fc1_plus10kb.bed -b peaks/chrom_curated/edu_B12_B10_peaks.bed -c > peaks/chrom_curated/count_edu_B12_B10_in_group_B10_up_fc1.bed")
system("bedtools intersect -a peaks/chrom_curated/groups_B10_down_fc1_plus10kb.bed -b peaks/chrom_curated/edu_B12_B10_peaks.bed -c > peaks/chrom_curated/count_edu_B12_B10_in_group_B10_down_fc1.bed")
system("bedtools intersect -a peaks/chrom_curated/groups_B10_up_fc1_plus10kb.bed -b peaks/chrom_curated/edu_B10_noB12_peaks.bed -c > peaks/chrom_curated/count_edu_B10_noB12_in_group_B10_up_fc1.bed")
system("bedtools intersect -a peaks/chrom_curated/groups_B10_down_fc1_plus10kb.bed -b peaks/chrom_curated/edu_B10_noB12_peaks.bed -c > peaks/chrom_curated/count_edu_B10_noB12_in_group_B10_down_fc1.bed")



# Leer el archivo con el conteo de picos por gen
counts.B12 <- read.table("peaks/intersect/count_edu_B12_in_genes.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-c(5,6)]
counts.B10 <- read.table("peaks/intersect/count_edu_B10_in_genes.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-c(5,6)]
counts.B2  <- read.table("peaks/intersect/count_edu_B2_in_genes.bed",  header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-c(5,6)]

counts.common.up <- read.table("peaks/chrom_curated/count_edu_B12_B10_in_group_B10_up_fc1.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-c(5,6)]
counts.common.down <- read.table("peaks/chrom_curated/count_edu_B12_B10_in_group_B10_down_fc1.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-c(5,6)]
counts.new.up <- read.table("peaks/chrom_curated/count_edu_B10_noB12_in_group_B10_up_fc1.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-c(5,6)]
counts.new.down <- read.table("peaks/chrom_curated/count_edu_B10_noB12_in_group_B10_down_fc1.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-c(5,6)]








# Asignar nombres de columnas
colnames(counts.B12) <- colnames(counts.B10) <- colnames(counts.B2) <- c("chr", "start", "end", "gene_id", "n_peaks")
colnames(counts.common.up) <- colnames(counts.common.down) <- colnames(counts.new.up) <- colnames(counts.new.down) <- c("chr", "start", "end", "gene_id", "n_peaks")


# Calcular longitud del gen
counts.B12$length <- counts.B12$end - counts.B12$start
counts.B10$length <- counts.B10$end - counts.B10$start
counts.B2$length  <- counts.B2$end  - counts.B2$start

counts.B12$sample <- "B12"
counts.B10$sample <- "B10"
counts.B2$sample  <- "B2"

counts.common.up$Peaks <- "Common"
counts.common.down$Peaks <- "Common"
counts.new.up$Peaks  <- "New"
counts.new.down$Peaks  <- "New"

# ----------- HISTOGRAM ALL GENES -----

counts <- rbind(counts.B12,counts.B10, counts.B2)
counts$n_peaks <- as.factor(counts$n_peaks)

# Crear tabla de frecuencia
library(dplyr)
df_freq <- counts %>%
  group_by(sample, n_peaks) %>%
  summarise(count = n()) %>%
  ungroup()

# Asegurar que n_peaks esté ordenado correctamente
df_freq$sample <- factor(df_freq$sample, levels = c("B12", "B10","B2"))

pdf("plots/boxplots/gene_clusters_length/histogram_gene_length_edu_mut.pdf")
ggplot(df_freq, aes(x = n_peaks, y = count, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("B12" = "steelblue", "B10" = "darkorange", "B2" = "seagreen")) +
  labs(title = "Distribución de número de picos por gen",
       x = "Número de picos por gen",
       y = "Número de genes") +
  theme_minimal()
dev.off()


# ----------- HISTOGRAM GROUPS -----

counts <- rbind(counts.common.down, counts.new.down)
counts$n_peaks <- as.factor(counts$n_peaks)

# Crear tabla de frecuencia
df_freq <- counts %>%
  group_by(Peaks, n_peaks) %>%
  summarise(count = n()) %>%
  ungroup()

# Asegurar que n_peaks esté ordenado correctamente
df_freq$Peaks <- factor(df_freq$Peaks, levels = c("Common", "New"))

pdf("plots/histogram/common_peaks_in_group_B10_down_fc1.pdf")
ggplot(df_freq, aes(x = n_peaks, y = count, fill = Peaks)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Common" = "#a5272b", "New" = "#d19883")) +
  labs(title = "Down group",
       x = "Number of peaks within gene",
       y = "Number of genes") +
  theme_minimal()
dev.off()

counts <- rbind(counts.common.up, counts.new.up)
counts$n_peaks <- as.factor(counts$n_peaks)

# Crear tabla de frecuencia
df_freq <- counts %>%
  group_by(Peaks, n_peaks) %>%
  summarise(count = n()) %>%
  ungroup()

# Asegurar que n_peaks esté ordenado correctamente
df_freq$Peaks <- factor(df_freq$Peaks, levels = c("Common", "New"))

pdf("plots/histogram/common_peaks_in_group_B10_up_fc1.pdf")
ggplot(df_freq, aes(x = n_peaks, y = count, fill = Peaks)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Common" = "#a5272b", "New" = "#d19883")) +
  labs(title = "Up group",
       x = "Number of peaks within gene",
       y = "Number of genes") +
  theme_minimal()
dev.off()



# ----------- VIOLINS --------

counts.B12$peak_group <- counts.B12$n_peaks
counts.B12$peak_group <- cut(counts.B12$n_peaks,
                             breaks = c(-1, 0, 3,5, Inf),
                             labels = c("0", "1–2", "3–5",">5"))
pdf('plots/boxplots/gene_clusters_length/edu_peaks_B12_genes_length.pdf')
ggplot(counts.B12, aes(x = peak_group, y = length)) +
  geom_violin(fill = "orchid", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  labs(title = "Distribución de longitud génica por grupo de picos en B12",
       x = "Grupo de número de picos",
       y = "Longitud del gen (bp)") +
  theme_minimal()
dev.off()

counts.B10$peak_group <- cut(counts.B10$n_peaks,
                             breaks = c(-1, 0, 5, 10, 15, 20, Inf),
                             labels = c("0", "1-5", "6-10", "11-15", "16-20",">20"))
pdf('plots/boxplots/gene_clusters_length/edu_peaks_B10_genes_length.pdf')
ggplot(counts.B10, aes(x = peak_group, y = length)) +
  geom_violin(fill = "orchid", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  labs(title = "Distribución de longitud génica por grupo de picos en B10",
       x = "Grupo de número de picos",
       y = "Longitud del gen (bp)") +
  theme_minimal()
dev.off()

counts.B2$peak_group <- cut(counts.B2$n_peaks,
                            breaks = c(-1, 0, 5, Inf),
                            labels = c("0", "1-5", ">5"))
pdf('plots/boxplots/gene_clusters_length/edu_peaks_B2_genes_length.pdf')
ggplot(counts.B2, aes(x = peak_group, y = length)) +
  geom_violin(fill = "orchid", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  labs(title = "Distribución de longitud génica por grupo de picos en B2",
       x = "Grupo de número de picos",
       y = "Longitud del gen (bp)") +
  theme_minimal()
dev.off()






############################## PEAKS DISTRIBUTION BY LENGTH ##############################



# EACH CLONE PEAKS

# Calcular longitud del gen
B12.peaks$width <- width(B12.peaks)
B10.peaks$width <- width(B10.peaks)
B2.peaks$width  <- width(B2.peaks)

summary(B12.peaks$width)
summary(B10.peaks$width)
summary(B2.peaks$width)


# Combinar todos los data.frames en uno solo
all_peaks <- rbind(data.frame(length = width(B12.peaks), sample = "B12"),
                   data.frame(length = width(B10.peaks), sample = "B10"),
                   data.frame(length = width(B2.peaks ), sample = "B2" ))
all_peaks$sample <- factor(all_peaks$sample, levels = c("B12", "B10", "B2"))


bins <- c(0, 1000, 5000, 10000, 20000, 50000, 100000, Inf)

# Asigna cada pico a su intervalo
all_peaks$length_bin <- cut(all_peaks$length,
                            breaks = bins,
                            labels = c("<1kb", "1–5kb", "5–10kb", "10–20kb", "20–50kb", "50–100kb", ">100kb"),
                            include.lowest = TRUE, right = FALSE)

# Verifica que la clasificación tenga sentido
table(all_peaks$length_bin, all_peaks$sample)

summary_df <- all_peaks %>%
  group_by(sample, length_bin) %>%
  summarise(count = n(), .groups = "drop")
total_peaks <- summary_df %>%
  group_by(sample) %>%
  summarise(total = sum(count))

# Unir los totales al dataframe resumen para calcular el porcentaje
summary_df <- summary_df %>%
  left_join(total_peaks, by = "sample") %>%
  mutate(percentage = (count / total) * 100)

# 🔹 Gráfico de barras con porcentajes
ggplot(summary_df, aes(x = length_bin, y = percentage, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  scale_fill_manual(values = c("B12"="#a5272b", "B10"="#d19883", "B2"="#6da1c2")) +
  labs(x = "Longitud de picos (bp)", y = "Porcentaje de picos (%)", 
       title = "Distribución de Longitud de Picos de EdU") +
  theme_minimal(base_size = 13) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# BOXPLOTS


# Combinar todos los data.frames en uno solo
all_peaks <- rbind(data.frame(length = width(B12.peaks      ), sample = "B12"      ),
                   data.frame(length = width(B10.peaks      ), sample = "B10"      ),
                   data.frame(length = width(B2.peaks       ), sample = "B2"       ),
                   data.frame(length = width(B12.B10.peaks  ), sample = "B12.B10"  ),
                   data.frame(length = width(B12.B2.peaks   ), sample = "B12.B2"   ),
                   data.frame(length = width(B10.B12.peaks  ), sample = "B10.B12"  ),
                   data.frame(length = width(B2.B12.peaks   ), sample = "B2.B12"   ),
                   data.frame(length = width(B12.noB10.peaks), sample = "B12.noB10"),
                   data.frame(length = width(B12.noB2.peaks ), sample = "B12.noB2" ),
                   data.frame(length = width(B10.noB12.peaks), sample = "B10.noB12"),
                   data.frame(length = width(B2.noB12.peaks ), sample = "B2.noB12" ))
all_peaks$sample <- factor(all_peaks$sample,
                           levels = c("B12", "B12.B10", "B12.B2", "B12.noB10", "B12.noB2",
                                      "B10", "B10.B12", "B10.noB12",
                                      "B2", "B2.B12", "B2.noB12"))



boxplot(all_peaks$length ~ all_peaks$sample, outline = F)

ggplot(all_peaks, aes(x = sample, y = length, fill = sample)) +
  geom_boxplot(outlier.shape = NA,  # Eliminar los outliers (puntos)
               color = "black",    # Color de las líneas de las cajas
               width = 0.6,        # Ajustar el ancho de las cajas
               alpha = 0.8) +      # Transparencia para mejor visualización
  scale_fill_manual(values = c("B12"= "#a5272b", "B12.B10"= "#a5272b", "B12.B2"= "#a5272b", "B12.noB10"= "#a5272b", "B12.noB2"= "#a5272b",
                               "B10"= "#d19883", "B10.B12"= "#d19883", "B10.noB12"= "#d19883",
                               "B2"= "#6da1c2", "B2.B12"= "#6da1c2", "B2.noB12"= "#6da1c2")) +  # Colores personalizados para cada muestra
  labs(x = "Muestras", y = "Longitud de Picos (bp)", 
       title = "Distribución de Longitud de Picos de EdU por Muestra") +
  theme_minimal(base_size = 14) +  # Estilo minimalista con un tamaño de base mayor
  theme(
    legend.title = element_blank(),  # Eliminar título de leyenda
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotar etiquetas del eje X
    axis.text.y = element_text(size = 12),  # Tamaño de las etiquetas del eje Y
    axis.title = element_text(size = 14),  # Tamaño de los títulos de los ejes
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Título del gráfico en negrita y centrado
    panel.grid.major = element_blank(),  # Eliminar la cuadrícula mayor
    panel.grid.minor = element_blank(),  # Eliminar la cuadrícula menor
    panel.border = element_blank()  # Eliminar borde del panel
  )+
  coord_cartesian(ylim = c(0, 80000))  # Ajustar el eje Y para que vaya de 0 a 150000




############################## CUSTOM ANNOTATION ##############################

# Data from /media/gmzlab/4TB/projects/h3k37/results/commands.sh
data <- data.frame(group=c("total", "genes", "genes05", "eugenes","eugenes05"),
                   common=c(5934, 4342, 3126, 3613, 2352),
                   new=c(3647, 2321, 2128, 1540, 1326))

data[-1] <-t(100*t(data[-1])/as.numeric(data[1,2:3]))


data_long <- data %>%
  pivot_longer(
    cols = c("common", "new"),
    names_to = c("type"),
    names_pattern = "(.*)",
    values_to = "value")%>%
  mutate(type = factor(type, levels = c("common", "new")))%>%
  mutate(group = factor(group, levels = c("total", "genes", "genes05", "eugenes","eugenes05")))
    


# Gráfico
pdf("plots/barplots/barplot_mass_spec_h3.pdf", height = 5, width = 25)
ggplot(data_long, aes(x = group, y = value, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  
  scale_fill_manual(values = c(
    "common" = "#a5272b", 
    "new" = "#d19883"    )) +
  labs(x="",y = "% peaks",
       fill = "type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

