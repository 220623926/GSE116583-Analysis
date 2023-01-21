# Installing packages
install.packages('BiocManager')
library(BiocManager)
install(c('tximport', 'DESeq2', 'biomaRt', 'pheatmap'))
library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)


# Download and read count data

sample_table = read_csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')
files = pull(sample_table, Run)
files = paste0('counts/', files, '/quant.sf')
names(files) = pull(sample_table, Run)
gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
txi = tximport(files, 
                 type='salmon',
                 tx2gene=gene_map,
                 ignoreTxVersion=TRUE)


# Applying normalisation & stats (Ward's binomial test)

dds = DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ Group)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
# the above 3 commands can be replaced by dds = DESeq(dds)
plotDispEsts(dds)



# Regularised log transformation - DO NOT use for DGE analysis
rld = rlog(dds)


# Principal component analysis (PCA) - Customised in gglot
pcaData <- plotPCA(rld, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c("blue", "orange","grey")) +
  theme_bw() +
  coord_fixed()


# Sample distance analysis - Euclidean distance
sample_distance = dist(t(assay(rld)), method='euclidian')
sample_distance_matrix = as.matrix(sample_distance)
heatmap_annotation = data.frame(group=colData(dds)[,c('Group')], row.names=rownames(colData(dds)))
ann_colors<-list(group=c(Naive="grey", Allo2h="orange", Allo24h="blue"))
pheatmap(sample_distance_matrix,
           clustering_distance_rows=sample_distance,
           clustering_distance_cols=sample_distance,
           width=7.3, height=5, color = brewer.pal(n = 9, name = 'BuPu'),
           annotation_colors = ann_colors,
           annotation_col = heatmap_annotation, filename = "Euclidian distance.png")


# DE Genes results

# Allo24h vs Naive
results_table = results(dds, contrast= c('Group', 'Allo24h', 'Naive'))
summary(results_table)
# Filter out NAs in rows:
results_tibble = as_tibble(results_table, rownames='ensembl_gene_id')
filtered_results = filter(results_tibble, complete.cases(results_tibble))

# Volcano plot (double-filtering)
filtered_results = mutate(filtered_results, logPVal = -log10(padj))

filtered_results = mutate(filtered_results, significant=ifelse(padj<0.05&abs(log2FoldChange)>1,
                                             'both', 
                                             ifelse(padj<0.05&abs(log2FoldChange)<1,
                                                    'pval only', 
                                                    'neither')))



ggplot(filtered_results, aes(x=log2FoldChange, y=logPVal)) +
  xlab(~log[2]~'FoldChange') +
  ylab(~-log[10]~'(Padj)') +
  xlim(-20,20) +
  geom_point(aes(colour=significant)) + 
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  geom_vline(xintercept=c(-1,1), linetype=2) +
  theme_bw()

ggsave("Allo24h vs Naive.png")



# Allo2h vs Naive
results_table = results(dds, contrast= c('Group', 'Allo2h', 'Naive'))
summary(results_table)
# Filter out NAs in rows:
results_tibble = as_tibble(results_table, rownames='ensembl_gene_id')
filtered_results = filter(results_tibble, complete.cases(results_tibble))
# Volcano plot (double-filtering)
filtered_results = mutate(filtered_results, logPVal = -log10(padj))

filtered_results = mutate(filtered_results, significant=ifelse(padj<0.05&abs(log2FoldChange)>1,
                                                               'both', 
                                                               ifelse(padj<0.05&abs(log2FoldChange)<1,
                                                                      'pval only', 
                                                                      'neither')))


ggplot(filtered_results, aes(x=log2FoldChange, y=logPVal)) +
  xlab(~log[2]~'FoldChange') +
  ylab(~-log[10]~'(Padj)') +
  xlim(-20,20) +
  geom_point(aes(colour=significant)) + 
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  geom_vline(xintercept=c(-1,1), linetype=2) +
  theme_bw()


ggsave("Allo2h vs Naive.png")



# Allo24h vs Allo2h
results_table = results(dds, contrast= c('Group', 'Allo24h', 'Allo2h'))
summary(results_table)
# Filter out NAs in rows:
results_tibble = as_tibble(results_table, rownames='ensembl_gene_id')
filtered_results = filter(results_tibble, complete.cases(results_tibble))
# Volcano plot (double-filtering)
filtered_results = mutate(filtered_results, logPVal = -log10(padj))


filtered_results = mutate(filtered_results, significant=ifelse(padj<0.05&abs(log2FoldChange)>1,
                                                               'both', 
                                                               ifelse(padj<0.05&abs(log2FoldChange)<1,
                                                                      'pval only', 
                                                                      'neither')))

ggplot(filtered_results, aes(x=log2FoldChange, y=logPVal)) +
  xlab(~log[2]~'FoldChange') +
  ylab(~-log[10]~'(Padj)') +
  xlim(-20,20) +
  geom_point(aes(colour=significant)) + 
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  geom_vline(xintercept=c(-1,1), linetype=2) +
  theme_bw()



ggsave("Allo24h vs Allo2h.png")


# Annotation of results table
ensembl108 = useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)
annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                  'start_position', 'end_position', 
                                  'strand', 'gene_biotype', 'external_gene_name', 
                                  'description'), 
                     filters = 'ensembl_gene_id', values = filtered_results$ensembl_gene_id, 
                     mart = ensembl108)
annot_results = left_join(filtered_results, annotation)
annot_results = arrange(annot_results, padj)
View(head(annot_results, 10))
degs = filter(annot_results, abs(log2FoldChange) > 1 & padj < 0.05)

degs_table = head(degs, 10)

write.csv(annot_results, file = "All DEGs - Allo24h vs Allo2h.csv")
write.csv(degs_table, file = "Top 10 DEGs - Allo24h vs Allo2h.csv")




