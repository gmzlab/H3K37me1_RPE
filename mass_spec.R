
############################## LOADING DATA AND FUNCTIONS ##############################

# Setting working directory and loading all necessary packages and data.

setwd('/media/gmzlab/4TB/projects/h3k37/results/')
#suppressMessages(library(tidyverse))
#suppressMessages(library(GenomicRanges))
#suppressMessages(library(gUtils))
suppressMessages(library(tidyverse))
#suppressMessages(library(ggupset))
#suppressMessages(library(GenomicRanges))
#suppressMessages(library(ggpubr))
#suppressMessages(library(cowplot))
#suppressMessages(library(heatmaply))
suppressMessages(library(ggplot2))
suppressMessages(library(genefilter))
suppressMessages(library(openxlsx))
#install.packages("pheatmap")
#suppressMessages(library(pheatmap))



############################## LOADING TABLE ##############################

data <- na.omit(read.xlsx('Gonzalo - Canonical PTMs - April 2025.xlsx', 1)[-1,c(1,3,5,7,2,4,6)])
colnames(data) <- c("name", "WT 1", "WT 2", "WT 3", "H3.3K37R 1", "H3.3K37R 2", "H3.3K37R 3")
data[,-1] <-lapply(data[,-1], as.numeric)
rownames(data) <- data$name
data <- data[!(data$name %in% c("H4_24_35 unmod","H4_79_92 unmod","H2A_82_88 unmod")),]
data <-data[rowSums(data[,-1])>0,]
#data$FC <- rowMeans(log2(data[,2:4]+1)) -rowMeans(log2(data[,5:7]+1))
data$FC <- log2(rowMeans(data[,5:7])) -log2(rowMeans(data[,2:4]))
data$p.value <- rowttests(as.matrix(data[,2:7]), factor(c(rep("WT", 3), rep("H3.3K37R", 3))) )$p.value
data$type <- "ns"
fc.threshold <- 1
pval.threshold <- 0.05
data$type[data$FC >  fc.threshold & data$p.value < pval.threshold] <- "up"
data$type[data$FC < -fc.threshold & data$p.value < pval.threshold] <- "down"

#onlyzero <- rowSums(data[,-1]) == 0
#data[,-1] <- t(apply(data[,-1], 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
#data[onlyzero,-1] <- 0

data.h3 <- data[grep("H3_", data$name),]
data.h4 <- data[grep("H4_", data$name),]





pdf("plots/heatmaps/mass_spec_h3.pdf", height = 20, width = 5)
pheatmap(
  data.h3[,2:7], 
  scale = "row",  # Escalar por filas (esto normaliza cada fila)
  color = colorRampPalette(c("darkblue", "skyblue","white", "#faff66", "red"))(100),  # Definir la paleta de colores
  cluster_rows = F,  # Hacer clustering en las filas (modificaciones de histona)
  cluster_cols = F,  # Hacer clustering en las columnas (muestras)
  show_rownames = TRUE,  show_colnames = TRUE
)
dev.off()


pdf("plots/heatmaps/mass_spec_h4.pdf", height = 10, width = 5)
pheatmap(
  data.h4[,2:7], 
  scale = "row",  # Escalar por filas (esto normaliza cada fila)
  color = colorRampPalette(c("darkblue", "skyblue","white", "#faff66", "red"))(100),  # Definir la paleta de colores
  cluster_rows = F,  # Hacer clustering en las filas (modificaciones de histona)
  cluster_cols = F,  # Hacer clustering en las columnas (muestras)
  show_rownames = TRUE,  show_colnames = TRUE
)
dev.off()


data.h3$color <- "ns"
data.h3$color[grep("K36me1",data.h3$name)] <- "K36me1"
data.h3$color[grep("K36me2",data.h3$name)] <- "K36me2"
data.h3$color[grep("K36me3",data.h3$name)] <- "K36me3"

# Orden de aparición (primero colores, luego grises)
data.h3$color <- factor(data.h3$color,
                        levels = c("K36me1", "K36me2", "K36me3", "ns"))

volcol <- c('#09635f','navy','darkgreen','grey')
names(volcol) <- c("K36me1", "K36me2", "K36me3", "ns")


highlight <- c("H3_27_40 K36me1", "H3_27_40 K36me2", "H3_27_40 K36me3","H3_27_40 K27me1K36me1","H3_27_40 K27me1K36me2","H3_27_40 K27me1K36me3","H3_27_40 K27me1K36me1","H3_27_40 K27me2K36me2","H3_27_40 K27me3K36me1","H3_27_40 K27me3K36me2")





pdf("plots/volcanos/volcano_mass_spec_h3.pdf", height = 5, width = 10)
ggplot(data.h3, aes(x = FC, y = -log10(p.value))) +
  
  # 1. Puntos grises (ns)
  geom_point(data = subset(data.h3, color == "ns"),
             aes(color = color),
             size = 3) +
  
  # 2. Puntos coloreados (K36me1/2/3)
  geom_point(data = subset(data.h3, color != "ns"),
             aes(color = color),
             size = 3) +
  
  scale_colour_manual(
    values = volcol,
    name = ""      # ← título de la leyenda (puedes cambiarlo)
  ) +
  
  ylim(0, 4) +  
  xlim(-6, 6) + 
  xlab('Log2 (Fold Change)') + 
  ylab('- Log10 (p-value)') +
  
  theme(
    legend.position = "right",   # ← ahora mostramos la leyenda
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    axis.line.x.bottom = element_line(color = 'black'),
    axis.line.y.left   = element_line(color = 'black'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 19)
  ) +
  
  annotate("segment", x = -fc.threshold, xend = -fc.threshold,
           y = -log10(pval.threshold), yend = 4,
           linetype = "dashed", color = "black", linewidth = 0.5) +
  
  annotate("segment", x = fc.threshold, xend = fc.threshold,
           y = -log10(pval.threshold), yend = 4,
           linetype = "dashed", color = "black", linewidth = 0.5) + 
  
  annotate("segment", x = -fc.threshold, xend = -6,
           y = -log10(pval.threshold), yend = -log10(pval.threshold),
           linetype = "dashed", color = "black", linewidth = 0.5) +
  
  annotate("segment", x = fc.threshold, xend = 6,
           y = -log10(pval.threshold), yend = -log10(pval.threshold),
           linetype = "dashed", color = "black", linewidth = 0.5)
dev.off()



data.h4$color <- "ns"

pdf("plots/volcanos/volcano_mass_spec_h4.pdf", height = 5, width = 10)
ggplot(data.h4, aes(x = FC, y = -log10(p.value))) +
  
  # 1. Puntos grises (ns)
  geom_point(data = subset(data.h4, color == "ns"),
             aes(color = color),
             size = 3) +
  
  # 2. Puntos coloreados (K36me1/2/3)
  geom_point(data = subset(data.h4, color != "ns"),
             aes(color = color),
             size = 3) +
  
  scale_colour_manual(
    values = volcol,
    name = ""      # ← título de la leyenda (puedes cambiarlo)
  ) +
  
  ylim(0, 4) +  
  xlim(-6, 6) + 
  xlab('Log2 (Fold Change)') + 
  ylab('- Log10 (p-value)') +
  
  theme(
    legend.position = "right",   # ← ahora mostramos la leyenda
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    axis.line.x.bottom = element_line(color = 'black'),
    axis.line.y.left   = element_line(color = 'black'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 19)
  ) +
  
  annotate("segment", x = -fc.threshold, xend = -fc.threshold,
           y = -log10(pval.threshold), yend = 4,
           linetype = "dashed", color = "black", linewidth = 0.5) +
  
  annotate("segment", x = fc.threshold, xend = fc.threshold,
           y = -log10(pval.threshold), yend = 4,
           linetype = "dashed", color = "black", linewidth = 0.5) + 
  
  annotate("segment", x = -fc.threshold, xend = -6,
           y = -log10(pval.threshold), yend = -log10(pval.threshold),
           linetype = "dashed", color = "black", linewidth = 0.5) +
  
  annotate("segment", x = fc.threshold, xend = 6,
           y = -log10(pval.threshold), yend = -log10(pval.threshold),
           linetype = "dashed", color = "black", linewidth = 0.5)
dev.off()



############################ BARPLOT ################

# ------------- H3
data_long <- data.h3 %>%
  pivot_longer(
    cols = c(`WT 1`, `WT 2`, `WT 3`, `H3.3K37R 1`, `H3.3K37R 2`, `H3.3K37R 3`),
    names_to = c("Genotype", "replica"),
    names_pattern = "(.*) (.*)",
    values_to = "valor"
  ) %>%
  mutate(
    Genotype = ifelse(grepl("WT", Genotype), "WT", "H3.3K37R")
  )
# Calcular medias por modificación y Genotype
data_summary <- data_long %>%
  group_by(name, Genotype) %>%
  summarise(media = mean(valor), sd = sd(valor), .groups = "drop") %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "H3.3K37R")))


# Gráfico
pdf("plots/barplots/barplot_mass_spec_h3.pdf", height = 5, width = 25)
ggplot(data_summary, aes(x = name, y = media, fill = Genotype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_manual(values = c(
    "WT" = "#a5272b",       # azul
    "H3.3K37R" = "#d19883"  # rojo
  )) +
  labs(x="",y = "Ratio",
       fill = "Genotype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# ------------- H4
data_long <- data.h4 %>%
  pivot_longer(
    cols = c(`WT 1`, `WT 2`, `WT 3`, `H3.3K37R 1`, `H3.3K37R 2`, `H3.3K37R 3`),
    names_to = c("Genotype", "replica"),
    names_pattern = "(.*) (.*)",
    values_to = "valor"
  ) %>%
  mutate(
    Genotype = ifelse(grepl("WT", Genotype), "WT", "H3.3K37R")
  )
# Calcular medias por modificación y Genotype
data_summary <- data_long %>%
  group_by(name, Genotype) %>%
  summarise(media = mean(valor), sd = sd(valor), .groups = "drop") %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "H3.3K37R")))


# Gráfico
pdf("plots/barplots/barplot_mass_spec_h4.pdf", height = 5, width = 15)
ggplot(data_summary, aes(x = name, y = media, fill = Genotype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_fill_manual(values = c(
    "WT" = "#a5272b",       # azul
    "H3.3K37R" = "#d19883"  # rojo
  )) +
  labs(x="",y = "Ratio",
       fill = "Genotype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


