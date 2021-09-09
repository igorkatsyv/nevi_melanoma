#### Melanoma scRNA-seq data
rm(list=ls())
library(mygene)
library(dplyr)
library(ggplot2)
library(ggExtra)

setwd('/Users/ikatsyv/Dropbox/nevi_melanoma')

# load data
expr = read.delim('GSE72056_melanoma_single_cell_revised_v2.txt')

expr2 = data.frame(t(expr))
colnames(expr2) = expr2[1,]
expr2 = expr2[-1,]
write.table(expr2,'GSE72056_melanoma_scRNaseq_transpose.txt', quote = F, row.names = F, sep = '\t')

remove(expr2)

expr = read.delim('GSE72056_melanoma_scRNaseq_transpose.txt')
row.names(expr) = paste(expr$tumor, expr$malignant.1.no.2.yes.0.unresolved., row.names(expr), sep = '|')

exclude = c(58, 67, 72, 74, 75) # samples with no melanoma cells

expr_malignant = expr[expr$malignant.1.no.2.yes.0.unresolved. == 2 &
                        expr$non.malignant.cell.type..1.T.2.B.3.Macro.4.Endo..5.CAF.6.NK. == 0 &  
                        !(expr$tumor %in% exclude),-c(1:3)] # remove non-malignant cells and tumors with no melanoma cells

pca <- prcomp(expr_malignant)
plot(pca$x)
ggplot(data.frame(pca$x), aes(x = PC1, y=PC2)) + geom_point() + theme_minimal()


# update gene symbols
genes = names(expr_malignant)

genes_update = queryMany(genes, species = 'human', scopes = 'symbol', fields = 'symbol')
genes_update_known = data.frame(id = genes_update$query, symbol = genes_update$symbol)
genes_update_known$symbol = ifelse(is.na(genes_update_known$symbol), genes_update_known$id, genes_update_known$symbol)


names(expr_malignant) = genes_update_known$symbol[match(names(expr_malignant), genes_update_known$id)]

plot(expr_malignant$PRAME)

prame = data.frame(PRAME = expr_malignant$PRAME, index = seq(1,length(expr_malignant$PRAME),1))

p1 = ggplot(prame, aes(x = index, y = PRAME, color = PRAME)) +
  geom_point() +
  theme_minimal() + 
  geom_hline(yintercept = median(expr_malignant$PRAME) - 3*mad(expr_malignant$PRAME), linetype = 'dashed', color = 'orange', size = 1.4) +
  theme(legend.position = 'none')

pdf('prame_single_melanoma.pdf', height = 6, width = 10)
ggMarginal(p1, type = 'density' , fill = 'grey20', color = 'grey20', margins = 'y') 
dev.off()

# Separate PRAME-low cases
expr_malignant_pramelow = data.frame(t(expr_malignant[expr_malignant$PRAME < (median(expr_malignant$PRAME) - 3*mad(expr_malignant$PRAME)),]))

# expr_malignant_pramelow = data.frame(t(expr_malignant[expr_malignant$PRAME == 0,]))

expr_malignant_pramelow$median = apply(expr_malignant_pramelow, 1, function(x) median(x))
expr_malignant_pramelow = expr_malignant_pramelow[order(-expr_malignant_pramelow$median),]

expr_malignant_pramelow[expr_malignant_pramelow$median > 0,]

# expr_malignant_pramelow = expr_malignant_pramelow[expr_malignant_pramelow$median > 6,]

expr_malignant_pramelow$iqr1 = apply(expr_malignant_pramelow, 1, function(x) summary(x)[2])
expr_malignant_pramelow = expr_malignant_pramelow[expr_malignant_pramelow$iqr1 > 4,]

pdf('prame_low_genes.pdf', height = 12, width = 18)
ggheatmap(expr_malignant_pramelow[,-c(ncol(expr_malignant_pramelow)-1,ncol(expr_malignant_pramelow))],
          fontsize_col = 5, fontsize_row = 5, dendrogram = 'none')
dev.off()

pramelow_genes = data.frame(gene = row.names(expr_malignant_pramelow))

write.table(pramelow_genes,'pramelow_genes.txt', quote = F, row.names = F, sep = '\t')

# Subset PRAME-low tumors

pramelowtumors = gsub('X','',
                      unique(do.call(cbind, strsplit(
                        names(expr_malignant_pramelow[,-c(ncol(expr_malignant_pramelow)-1,ncol(expr_malignant_pramelow))]), 
                        '.', fixed = T))[1,]))

expr2 = data.frame(tumor = row.names(expr_malignant), PRAME = expr_malignant$PRAME)

expr2$Tumor = do.call(cbind, strsplit(expr2$tumor, '|', fixed = T))[1,]

expr2 = expr2[expr2$Tumor %in% pramelowtumors,]

expr2$Tumor = as.factor(expr2$Tumor)

pdf('pramelowtumors_pramexpr.pdf', height = 4, width = 8)
ggbetweenstats(expr2, x = Tumor, y = PRAME, plot.type = 'violin', 
               centrality.type = 'nonparametric' , 
               results.subtitle = F, pairwise.comparisons = F,
               palette = 'Set3') +geom_hline(yintercept = median(expr_malignant$PRAME) - 3*mad(expr_malignant$PRAME), linetype = 'dashed')
dev.off()


# # load DEGs
# melanoma_nevi_deg = read.delim('GSE112509_melanoma_nevi_deg.txt')
# 
# melanoma_nevi_deg_scfiltered = melanoma_nevi_deg[melanoma_nevi_deg$gene %in% row.names(expr_malignant_pramelow),]
# 
# write.table(melanoma_nevi_deg_scfiltered[melanoma_nevi_deg_scfiltered$logFC > 0,],'melanoma_nevi_deg_scfiltered.txt', quote = F, row.names = F, sep = '\t')

save.image('melanoma_scRNAseq.RData')