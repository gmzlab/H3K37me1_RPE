
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


############################## RARE CHROMOSOMES REMOVAL ##############################

peaks <- readPeakFile("peaks/h3k37me1_async/h3k37me1_async_peaks_fc3.narrowPeak")
peaks <- peaks[seqnames(peaks) %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")]
seqlevels(peaks) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
rtracklayer::export(peaks, "peaks/chrom_curated/h3k37me1_async_peaks_fc3.bed")




############################## PEAKS ANNOTATION ##############################

peak.file <- readPeakFile('peaks/chrom_curated/h3k37me1_async_peaks_fc3.bed')
pdf.file  <- 'plots/piecharts/piechart_peaks_h3k37me1_async_fc3_complex.pdf'

df <- as.data.frame(annotatePeak(peak.file, tssRegion=c(-1000, 1000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level = "gene", overlap = "all"))
#df <- as.data.frame(annotatePeak(peak.file, tssRegion=c(-1000, 1000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level = "gene", overlap = "all", genomicAnnotationPriority= c("5UTR", "3UTR", "Exon", "Intron", "Promoter", "Downstream", "Intergenic")))
df$region <- sub(" \\(.*", "", df$annotation)
df$group <- ifelse(df$region %in%  c("Intron", "Exon", "3' UTR", "5' UTR"), "Gene body",
            ifelse(df$region %in% c("Distal Intergenic", "Downstream"), "Intergenic", df$region))
df_plot <- df %>% dplyr::count(region) %>% mutate(porcentaje = round(100 * n / sum(n), 1),
           etiqueta = paste0(region, "\n", porcentaje, "%"))
orden_regiones <- c(
  setdiff(unique(df_plot$region),
          c("Distal Intergenic", "Promoter")),
  "Distal Intergenic",
  "Promoter"
)

df_plot$region <- factor(df_plot$region, levels = orden_regiones)
df_plot <- df_plot[order(df_plot$region), ]

# ---- Paleta de colores ----
colores <- c(
  "#2F8036",  "#429A48" ,"#57A85C",   "#6FBF73",    # regiones intermedias
  "#e39610",  # Distal Intergenic (azul)
  "#4f95ff"   # Promoter (rojo)
)

# ---- Plot ----
pdf(pdf.file)
ggplot(df_plot, aes(x = "", y = porcentaje, fill = region)) +
  geom_col(color = "white", width = 1) +
  coord_polar("y") +
  theme_void() +
  geom_text(aes(label = etiqueta),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = colores) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()
plotAnnoBar(annotatePeak(peak.file, tssRegion=c(-1000, 1000), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, level = "gene", overlap = "all"))

############################## GENES ##############################

system("bedtools intersect -u -a /home/gmzlab/Documents/ricardo/annotation/Homo_sapiens_ann_genes_112.bed -b peaks/chrom_curated/h3k37me1_async_peaks_fc3.bed -f 0.5 > peaks/chrom_curated/h3k37me1_async_genes_fc3.bed")

go <- clusterProfiler::enrichGO(gene = readPeakFile('peaks/chrom_curated/h3k37me1_async_genes_fc3.bed')$V4, OrgDb = org.Hs.eg.db, ont = 'ALL', 
      pvalueCutoff = 0.05, qvalueCutoff = 0.05, keyType = 'ENSEMBL')
go$Description



