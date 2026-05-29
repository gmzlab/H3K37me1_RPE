############### WD AND PACKAGES ####

setwd("/media/gmzlab/4TB/projects/h3k37/results/")

# Loading libraries

#install.packages("ggrepel")
suppressMessages(library(ggrepel))

suppressMessages(library(dplyr))
#BiocManager::install("ballgown")
suppressMessages(library(ballgown))
#BiocManager::install("clusterProfiler")
suppressMessages(library(clusterProfiler))
#install.packages("cowplot")
suppressMessages(library(cowplot))
#BiocManager::install("DESeq2")
suppressMessages(library(DESeq2))
#BiocManager::install("genefilter")
suppressMessages(library(genefilter))
#install.packages("factoextra")
suppressMessages(library(factoextra))
#install.packages("FactoMineR")
suppressMessages(library(FactoMineR))
#BiocManager::install("GenomicRanges")
suppressMessages(library(GenomicRanges))
#install.packages("ggplot2")
suppressMessages(library(ggplot2))
#install.packages("ggpubr")
suppressMessages(library(ggpubr))
#install.packages("htmltools")
suppressMessages(library(htmltools))
#BiocManager::install("limma")
suppressMessages(library(limma))
#install.packages("MetBrewer")
suppressMessages(library(MetBrewer))
#install.packages("openxlsx")
suppressMessages(library(openxlsx))
#BiocManager::install("org.Hs.eg.db")
suppressMessages(library(org.Hs.eg.db))
#BiocManager::install("org.Mm.eg.db")
suppressMessages(library(org.Mm.eg.db))
#BiocManager::install("org.Sc.sgd.db")
suppressMessages(library(org.Sc.sgd.db))
#install.packages("plotly")
suppressMessages(library(plotly))
#install.packages("tidyverse")
suppressMessages(library(tidyverse))
#BiocManager::install("fgsea")
suppressMessages(library(fgsea))
#install.packages(data.table)
suppressMessages(library(data.table))
#BiocManager::install('pathview')
suppressMessages(library("pathview"))



# Loading all annotated features.
#ann <- rtracklayer::import('/home/gmzlab/Documents/ricardo/annotation/Homo_sapiens.GRCh38.112.gtf')
# Getting all annotated genes by generating a GRanges object with the corresponding m columns
#write.table(as.matrix(mcols(ann[ann$type == 'gene'])[c('gene_id', 'gene_name')]), '/home/gmzlab/Documents/ricardo/annotation/geneid2name.Homo_sapiens.GRCh38.112.csv', sep = ',', quote = F, row.names = F)
geneid2name <- read.csv2('/home/gmzlab/Documents/ricardo/annotation/geneid2name.Homo_sapiens.GRCh38.112.csv', sep = ',')


experimental.design <- data.frame(sample=c('chrna_B2_1','chrna_B2_2','chrna_B10_1','chrna_B10_2','chrna_B12_1','chrna_B12_2'),condition=factor(c('mut2', 'mut2', 'mut1', 'mut1', 'wt', 'wt')))
experimental.design <- data.frame(sample=c('chrna_B2_1','chrna_B2_2','chrna_B12_1','chrna_B12_2'),condition=factor(c('mut', 'mut', 'wt', 'wt')))
experimental.design <- data.frame(sample=c('chrna_B10_1','chrna_B10_2','chrna_B12_1','chrna_B12_2'),condition=factor(c('mut', 'mut', 'wt', 'wt')))
experimental.design <- experimental.design[sort(experimental.design$sample,ind=T)$ix,]
number.samples <- nrow(experimental.design)
number.replicates <- table(experimental.design$condition)
sample.labels <- experimental.design$sample




############### DESEQ2 ####

control.condition <- 'wt'
experimental.condition <- 'mut2'


control.indeces <- vector(mode="numeric",length=number.replicates[[control.condition]])
experimental.indeces <- vector(mode="numeric",length= number.replicates[[experimental.condition]])
j <- 1
k <- 1
for(i in 1:nrow(experimental.design))
{
  if(experimental.design$condition[i] == control.condition)
  {
    control.indeces[j] <- i
    j <- j + 1
  } else if (experimental.design$condition[i] == experimental.condition)
  {
    experimental.indeces[k] <- i
    k <- k + 1
  }
}

contrast.experimental.design <- experimental.design[c(control.indeces,experimental.indeces),]
contrast.experimental.design$condition <- as.factor(contrast.experimental.design$condition)

## Loading data with raw counts
#gene.count.matrix <- as.matrix(read.csv(file="quantification_stringtie/gene_count_matrix_edu.csv", sep=",", header = T, comment.char = "#", row.names = 1)[,-c(1:5)])
#colnames(gene.count.matrix)<-gsub(".bam", "", colnames(gene.count.matrix))
#gene.count.matrix<- gene.count.matrix[,contrast.experimental.design$sample]

# desde strintie
#gene.count.matrix <- as.matrix(read.table(file = 'quantification_stringtie/gene_count_matrix_feature.csv', header = T, sep = ',', row.names = 1))[,contrast.experimental.design$sample] # Quitamos muestras
#row.names(gene.count.matrix) <- gsub('\\|.*','',row.names(gene.count.matrix))
#gene.count.matrix$names <- rownames(gene.count.matrix)
#gene.count.matrix <- as.matrix(gene.count.matrix)[-nrow(gene.count.matrix),]

# desde FEATURECOUNTS
gene.count.matrix <- as.matrix(read.table(file = 'quantification_stringtie/gene_count_matrix_feature.csv', header = T, sep = '\t', row.names = 1)[6:11])
colnames(gene.count.matrix) <- c("chrna_B10_1", "chrna_B10_2", "chrna_B12_1", "chrna_B12_2", "chrna_B2_1", "chrna_B2_2")
gene.count.matrix <- gene.count.matrix[,contrast.experimental.design$sample]



gene.count.matrix["ENSG00000181222",]
summary(gene.count.matrix)

gene.count.matrix <- gene.count.matrix[rowSums(gene.count.matrix>0) == ncol(gene.count.matrix),]
gene.count.df <- left_join(data.frame(gene_id = rownames(gene.count.matrix), gene.count.matrix), geneid2name, 'gene_id')[c('gene_id', 'gene_name', colnames(gene.count.matrix))]

dds <- DESeqDataSetFromMatrix(countData=gene.count.matrix, colData=contrast.experimental.design, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = control.condition)
dds <- DESeq(dds)

fc.threshold <- log2(2)
qval.threshold <- 0.05

de.results <- as.data.frame(results(dds)) 
de.results$gene_type <- 'ns'
de.results[which(de.results$log2FoldChange >  fc.threshold & de.results$padj < qval.threshold),'gene_type'] <- 'activated'
de.results[which(de.results$log2FoldChange < -fc.threshold & de.results$padj < qval.threshold),'gene_type'] <- 'repressed'
de.results$gene_id <- rownames(de.results)
rownames(de.results) <- NULL
de.results <- left_join(de.results, geneid2name, 'gene_id')
de.results <- de.results[c('gene_id','gene_name','baseMean','log2FoldChange','pvalue','padj', 'stat','gene_type')]
de.results <- de.results[de.results$baseMean > 0,]
de.results <- de.results[order(de.results$gene_type, -de.results$log2FoldChange, decreasing = F),]
de.results[,"link"] <- paste0('= HYPERLINK("https://www.genecards.org/cgi-bin/carddisp.pl?gene=', de.results[,"gene_id"], '")')
de.results <- de.results[!is.na(de.results$padj),]

table(de.results$gene_type)


#de.results <- read.xlsx('rnaseq_results.xlsx', 1)

wb <- createWorkbook()


addWorksheet(wb, paste0('DESEQ2_',experimental.condition, '_', control.condition))
writeData(wb, paste0('DESEQ2_',experimental.condition, '_', control.condition), de.results)
writeFormula(wb, paste0('DESEQ2_',experimental.condition, '_', control.condition), de.results[,'link'], startCol = ncol(de.results), startRow = 2)

addWorksheet(wb, paste0('activated_', experimental.condition, '_', control.condition))
writeData(wb, paste0('activated_', experimental.condition, '_', control.condition), de.results[de.results$gene_type=='activated',c('gene_id','gene_name')])

addWorksheet(wb, paste0('repressed_', experimental.condition, '_', control.condition))
writeData(wb, paste0('repressed_', experimental.condition, '_', control.condition), de.results[de.results$gene_type=='repressed',c('gene_id','gene_name')])

write.table(file="../../rnaseq_results_122024/repressed_feature.txt",de.results[de.results$gene_type=='repressed',c('gene_id')], quote = FALSE,row.names = F, col.names = F)





highlight.genes <- c("RNU6-264P")

pdf(paste0('plots/volcanos/volcano_', experimental.condition, '_', control.condition, '.pdf'))
volcol <- c('#e51a1f', '#204b9b','grey')
names(volcol) <- c("activated","repressed","ns")
ggplot(de.results, aes(x=log2FoldChange , y=-log10(padj), color=gene_type)) +
  geom_point(size = 1) +
  geom_point(data = subset(de.results, gene_name %in% highlight.genes), color = "black",size = 2) +
  geom_text_repel(data = subset(de.results, gene_name %in% highlight.genes),
                  aes(label = gene_name),
                  size = 4, max.overlaps = 10, colour = 'black') +
  scale_colour_manual(values = volcol) + 
  #ylim(0,10) +  
  #xlim(-5.5,5.5)+ 
  xlab('Log2 (Fold Change)') + ylab('- Log10 (q-Value)') +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line.x.bottom = element_line(color = 'black'),
        axis.line.y.left   = element_line(color = 'black'),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 19))+
  annotate("segment", x = -fc.threshold, xend = -fc.threshold, y = -log10(qval.threshold), yend=max(-log10(de.results$padj)), # 
           linetype = "dashed", color = "black", linewidth = 0.5) +
  annotate("segment", x = fc.threshold, xend = fc.threshold, y = -log10(qval.threshold), yend=max(-log10(de.results$padj)), # 
           linetype = "dashed", color = "black", linewidth=0.5) + 
  annotate("segment",x = -fc.threshold, xend = min(de.results$log2FoldChange), y = -log10(qval.threshold), yend=-log10(qval.threshold),
           linetype = "dashed", color = "black", linewidth=0.5) +
  annotate("segment",x = fc.threshold, xend = max(de.results$log2FoldChange), y = -log10(qval.threshold), yend=-log10(qval.threshold),
           linetype = "dashed", color = "black", linewidth=0.5)
dev.off()








############### GENE ONTOLOGY ############

go.activated <- enrichGO(gene = de.results$gene_id[de.results$gene_type=='activated'], OrgDb = org.Hs.eg.db, ont = 'BP', 
                         pvalueCutoff = 0.05, qvalueCutoff = 0.05, keyType = 'ENSEMBL')

go.activated$Description


addWorksheet(wb, paste0('GOactv_', experimental.condition, '_', control.condition))
writeData(wb, paste0('GOactv_', experimental.condition, '_', control.condition), as.data.frame(go.activated))

addWorksheet(wb, paste0('genes_GOactv_', experimental.condition, '_', control.condition))
writeData(wb, paste0('genes_GOactv_', experimental.condition, '_', control.condition), left_join(data.frame(gene_id = unique(unlist(strsplit(as.data.frame(go.activated)$geneID, split = '/')))),geneid2name,'gene_id'))


pdf(file='plots/dotplot_go_activated.pdf')
clusterProfiler::dotplot(go.activated, showCategory=c(
  "pyridine-containing compound catabolic process",
  "NADP metabolic process",
  "glutamine family amino acid catabolic process"
), x='GeneRatio', orderBy='x')
dev.off()


go.repressed <- enrichGO(gene = de.results$gene_id[de.results$gene_type=='repressed'], OrgDb = org.Hs.eg.db, ont = 'BP', 
                         pvalueCutoff = 0.05, qvalueCutoff = 0.05, keyType = 'ENSEMBL')

go.repressed$Description


addWorksheet(wb, paste0('GOrepr_', experimental.condition, '_', control.condition))
writeData(wb, paste0('GOrepr_', experimental.condition, '_', control.condition), as.data.frame(go.repressed))

addWorksheet(wb, paste0('genes_GOrepr_', experimental.condition, '_', control.condition))
writeData(wb, paste0('genes_GOrepr_', experimental.condition, '_', control.condition), left_join(data.frame(gene_id = unique(unlist(strsplit(as.data.frame(go.repressed)$geneID, split = '/')))),geneid2name,'gene_id'))


pdf(file='plots/dotplot_go_repressed.pdf')
clusterProfiler::dotplot(go.repressed, showCategory=10, x='GeneRatio', orderBy='x')
dev.off()



saveWorkbook(wb, paste0('rnaseq_results_', experimental.condition,'_',control.condition,'.xlsx'), overwrite = TRUE)


go.genes <- as.data.frame(go.activated)[, c('Description', 'geneID')]
go.genes <- as.data.frame(go.repressed)[, c('Description', 'geneID')]

wb.go <- createWorkbook()
term <- 'serine family amino acid biosynthetic process'
addWorksheet(wb.go, gsub( ':', '_', rownames(go.genes[go.genes$Description == term,])))
writeData(wb.go, gsub( ':', '_', rownames(go.genes[go.genes$Description == term,])),
          left_join(data.frame(gene_id=strsplit(go.genes[go.genes$Description == term, 'geneID'], '/')[[1]]), geneid2name, 'gene_id'))
saveWorkbook(wb.go, paste0('tables/', gsub(' ', '_', gsub(', ', '_', term)), '.xlsx'), overwrite = TRUE)





############### GSEA ####

go_pathways <- gmtPathways("/home/gmzlab/Documents/ricardo/annotation/c5.all.v2024.1.Hs.symbols.gmt")
gene.stat <- de.results %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% group_by(gene_name) %>% summarize(stat=mean(stat))

fgseaRes <- fgsea(pathways=go_pathways, stats=deframe(gene.stat))

# Show in a nice table:
go_gseaTidy <- as_tibble(fgseaRes) %>% arrange(desc(NES))


#wb <- createWorkbook()

addWorksheet(wb, paste0('GSEA_', experimental.condition, '_', control.condition))
writeData(wb, paste0('GSEA_', experimental.condition, '_', control.condition), go_gseaTidy)

saveWorkbook(wb, 'tables/gsea.xlsx', overwrite = TRUE)



plotEnrichment(go_pathways[["GOBP_RNA_SPLICING"]],
               deframe(gene.stat)) + labs(title="RNA Splicing")



filtered_terms <-c("GOBP_ALPHA_AMINO_ACID_METABOLIC_PROCESS",
                   "GOBP_MRNA_METABOLIC_PROCESS",
                   "GOBP_RNA_PROCESSING",
                   "GOBP_MYELOID_LEUKOCYTE_DIFFERENTIATION",
                   "GOBP_MYELOID_PROGENITOR_CELL_DIFFERENTIATION",
                   "GOBP_RNA_SPLICING",
                   "GOBP_TRNA_METABOLIC_PROCESS",
                   'GOBP_REGULATION_OF_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS')


go_filtered <- go_gseaTidy[go_gseaTidy$pathway %in% filtered_terms,] 
ggplot(go_filtered, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO pathways NES from GSEA") 



################# VENN DIAGRAM ####

rip <- readxl::read_excel("../ripseq/tables/rip_results.xlsx", sheet = 2)$gene_id
metabolic <- readxl::read_excel("tables/rnaseq_results.xlsx", sheet = 7)$gene_id
activated <- readxl::read_excel("tables/rnaseq_results.xlsx", sheet = 2)$gene_id
repressed <- readxl::read_excel("tables/rnaseq_results.xlsx", sheet = 3)$gene_id
activated.mouse <- read.table(file="tables/activated_genes_mouse_converted.txt", header=FALSE)$V1
repressed.mouse <- read.table(file="tables/repressed_genes_mouse_converted.txt", header=FALSE)$V1
majiq <- unique(gsub('gene:', '', read_delim("majiq/modulize_output/summary.tsv", delim = "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)$gene_id))

BioVenn::draw.venn(list_x = majiq, 
                   list_y = majiq.jeff,
                   list_z = NULL,
                   title = '', 
                   t_s = 50,
                   t_f= 'verdana',
                   subtitle = '', 
                   xtitle = 'Human',
                   xt_s = 3, xt_f = 'verdana',
                   ytitle = 'Mouse', 
                   yt_s = 3, yt_f = 'verdana',
                   ztitle = 'Repressed', 
                   zt_s = 35, zt_f = 'verdana',
                   output = 'tif', 
                   filename = 'plots/venn_majiq.tiff', 
                   width = 1000, 
                   height = 1000, 
                   nrtype = NULL,#nr_s = 10, nr_f = 'verdana', 
                   x_c = '#42acd7',
                   y_c='#e51a1f',
                   z_c='#204b9b'
)


'#e3ded3'
'#b79d71'
'#42acd7'
'#5f5f5f'
'#2d59a2'
'#204b9b'
'#e51a1f'


###### TEST EXACTO DE FISHER #######
tabla <- matrix(c(2, 199, 224, 14059), nrow = 2,
                dimnames = list(Grupo = c("Activated","No activated" ),
                                Resultado = c( "Splicing","No splicing")))

tabla <- matrix(c(4, 197, 309, 13974 ), nrow = 2,
                dimnames = list(Grupo = c("RIP","TOTAL" ),
                                Resultado = c( "Splicing","No splicing")))


fisher.test(tabla)


