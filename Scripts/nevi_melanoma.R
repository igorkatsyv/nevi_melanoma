#### Identify PRAME DiffCorr genes in melanocytic nevi vs. melanoma
rm(list=ls())

library(magrittr)
library(genefilter)
library(DGCA)
library(gplots)
library(mygene)
library(WGCNA)
library(limma)
library(data.table)
library(ggplot2)
library(EnhancedVolcano)
library(enrichR)

setwd('/Users/ikatsyv/Dropbox/nevi_melanoma')

# load expression data
# load expression data
expr = read.delim('GSE112509_DESeq2_normalized_counts.txt')
row.names(expr) = expr[,1] ; expr = expr[,-1]
row.names(expr) = do.call(cbind,strsplit(as.character(row.names(expr)), '.',fixed=T))[1,]

## Preprocess expr
expr = expr[rowSums(expr==0)<=5,] # Remove genes without non-zero counts in at least 5 samples
expr = log2(expr+1)

# remove genes with low variance
options(stringsAsFactors = FALSE)

datExpr0 = as.data.frame(t(expr));

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf('GSE112509_nevi_melanoma_sampleclustering.pdf',width=15, height = 8.8)
plot(sampleTree, main = "Sample heirarchical clustering by average Euclidian distance", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# Plot a line to show the cut
pdf('GSE112509_nevi_melanoma_sampleclustering_treecut.pdf',width=15, height = 8.8)
plot(sampleTree, main = "Sample heirarchical clustering by average Euclidian distance", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 288, col = "red");
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 288, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

expr = data.frame(t(datExpr))

# assign gene symbols to geneIDs
genes = row.names(expr)

genes_update = queryMany(genes, species = 'human', scopes = 'ensembl.gene', fields = 'symbol')
genes_update_known = data.frame(id = genes_update$query, symbol = genes_update$symbol)

genes_update_known = genes_update_known[complete.cases(genes_update_known),]

expr = expr[row.names(expr) %in% genes_update_known$id,]


expr$updated = genes_update_known$symbol[match(row.names(expr), genes_update_known$id)]
row.names(expr) = paste(row.names(expr),expr$updated, sep= '|')

expr = expr[,-76]

write.table(expr, 'GSE112509_melanoma_nevi_RNAseq_preprocessed.txt', quote = F, row.names = F, sep = '\t')


# Call DEGs
design = data.frame(samples = colnames(expr))
design$nevus = ifelse(grepl('_N',design$samples), 1, 0)
design$melanoma = ifelse(grepl('_M',design$samples), 1, 0)
rownames(design) = design[,1] ; design = design[,-1]

sml = c(rep("Nevus",nrow(design[design$nevus == 1,])), rep("Melanoma",nrow(design[design$melanoma == 1,])))

fl = as.factor(sml)
design = model.matrix(~ fl + 0, expr)
colnames(design) = levels(fl)
fit = lmFit(expr, design)
cont.matrix = makeContrasts(Melanoma-Nevus, levels=design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, adjust="fdr", sort.by="B", number=nrow(expr))
tT$gene = do.call(cbind,strsplit(row.names(tT),"|",fixed=T))[2,]
melanoma_nevus_deg = tT[,c(7,1:6)]
melanoma_nevus_deg$module = ifelse(melanoma_nevus_deg$logFC > 0, 'UP_melanoma', 'DN_melanoma')
melanoma_nevus_deg = melanoma_nevus_deg[2^(abs(melanoma_nevus_deg$logFC)) > 2 & melanoma_nevus_deg$adj.P.Val < 0.1,]
write.table(melanoma_nevus_deg,'GSE112509_melanoma_nevi_deg.txt', quote = F, row.names = F, sep = '\t')

expr_deg = expr[row.names(expr) %in% row.names(melanoma_nevus_deg),]

# load prame-low genes
pramelow = read.delim('pramelow_genes.txt')

ttpramelow = tT[tT$gene %in% pramelow$gene,]
ttpramelow = ttpramelow[ttpramelow$logFC > 0,]

pdf('pramelow_nm_deg.pdf', width = 16, height = 6.8)
EnhancedVolcano(ttpramelow, x = 'logFC', y = 'P.Value', lab = ttpramelow$gene, FCcutoff = 0, pCutoff = 1e-02)
dev.off()

ttpramelow$adj.P.Val = p.adjust(ttpramelow$P.Value, method = 'fdr')


ttpramelow = ttpramelow[ttpramelow$adj.P.Val < 0.01,]



# load human protein atlas data
hpa = read.delim('hpa_normal_tissue.tsv')

hpa_genes = unique(hpa$Gene)

hpa_genes_update = queryMany(hpa_genes, species = 'human', scopes = 'ensembl.gene', fields = 'symbol')

hpa_genes_update2 = data.frame(id = hpa_genes_update@listData$query,
                               symbol = hpa_genes_update@listData$symbol)


hpa$updated = hpa_genes_update2$symbol[match(hpa$Gene, hpa_genes_update2$id)]

hpa_backgroundhigh = hpa %>%
  filter(Tissue == 'skin 1') %>%
  filter(Cell.type == 'melanocytes' | Cell.type == 'keratinocytes') %>%
  filter(Level == 'Medium' | Level == 'High') %>% 
  filter(Reliability != 'Uncertain')


ttpramelow_hpafiltered = ttpramelow[!(ttpramelow$gene %in% hpa_backgroundhigh$updated),]

expr_candidates = expr[row.names(expr) %in% row.names(ttpramelow_hpafiltered),]
row.names(expr_candidates) = do.call(cbind, strsplit(row.names(expr_candidates), '|', fixed =T))[2,]

expr_candidates_class = data.frame(class = do.call(cbind,strsplit(names(expr_candidates), '_', fixed = T))[2,])
expr_candidates_class$color = ifelse(expr_candidates_class$class == 'N','deepskyblue','magenta')

clusters = kmeans(t(expr_candidates), 2)

pdf('pramelow_candidates.pdf', height = 4, width = 4)
heatmap(as.matrix(expr_candidates), cexRow = .5, ColSideColors = expr_candidates_class$color, labCol = '', distfun = function(x) dist(x, method = "manhattan"))
dev.off()

hpa_candidates = hpa[hpa$updated %in% ttpramelow_hpafiltered$gene &
                       hpa$Level == 'High' &
                       hpa$Reliability != 'Uncertain',]

# load HPA RNA data

hpa_rna = read.delim('hpa_rna_single_cell_type.tsv')

hpa_rna_genes = unique(hpa_rna$Gene)

hpa_rna_genes_updated = queryMany(hpa_rna_genes, species = 'human', scopes = 'ensembl.gene', fields = 'symbol')

hpa_rna_genes_updated2 = data.frame(id = hpa_rna_genes_updated$query, gene = hpa_rna_genes_updated$symbol)

hpa_rna$updated = hpa_rna_genes_updated2$gene[match(hpa_rna$Gene, hpa_rna_genes_updated2$id)]

hpa_rna_candidates = hpa_rna[hpa_rna$updated %in% ttpramelow_hpafiltered$gene,]



# perform functional annotations

pramelow_enrichr = enrichr(pramelow$gene, databases = c('MSigDB_Hallmark_2020', 'BioPlanet_2019','GeneSigDB','GO_Biological_Process_2021','KEGG_2021_Human','WikiPathway_2021_Human'))[[1]]
pramelow_enrichr$Term = factor(pramelow_enrichr$Term, levels = rev(pramelow_enrichr$Term))

# plotEnrich(pramelow_enrichr, title = 'Pathways enriched in PRAME-low melanocytes', order = 'P.Value', showTerms = 20)

pdf('pramelow_enrichr.pdf', height = 8, width = 8)
ggplot(pramelow_enrichr[1:20,], aes(x = Term, y = -log10(Adjusted.P.value), fill = Odds.Ratio)) +
  geom_bar(stat = 'identity') + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12)) +
  labs(fill = 'Odds Ratio') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'orange', size = 1)+
  coord_flip()
dev.off()

final_enrichr = enrichr(ttpramelow_hpafiltered$gene, databases = c('MSigDB_Hallmark_2020', 'BioPlanet_2019','GeneSigDB','GO_Biological_Process_2021','KEGG_2021_Human','WikiPathway_2021_Human'))[[1]]
final_enrichr$Term = factor(final_enrichr$Term, levels = rev(final_enrichr$Term))


# plotEnrich(final_enrichr, title = 'Pathways enriched in candidate genes', order = 'P.Value', showTerms = 20)

pdf('final_enrichr.pdf', height = 8, width = 8)
ggplot(final_enrichr[1:20,], aes(x = Term, y = -log10(Adjusted.P.value), fill = Odds.Ratio)) +
  geom_bar(stat = 'identity') + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12)) + 
  labs(fill = 'Odds Ratio') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'orange', size = 1)+
  coord_flip()
dev.off()

save.image('nevi_melanoma.RData')