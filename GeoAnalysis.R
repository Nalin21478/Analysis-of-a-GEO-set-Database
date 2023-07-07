library(GEOquery)
library(ggplot2)
library(tidyverse)
library(Biobase)
library(corrplot)
library(preprocessCore)

my_id <- "GSE16015"
gse <- getGEO(my_id)[[1]]
exprs_data <- exprs(gse)
dim(exprs_data)
summary(exprs(gse))
exprs_df <- as.data.frame(exprs_data)

cor_matrix <- cor(exprs_data)
corrplot(cor_matrix, method = "color")
boxplot(exprs(gse),outline=FALSE)

dim(exprs_data)
summary(exprs(gse))

exprs_data_norm <- normalize.quantiles(exprs_data)
exprs_df <- as.data.frame(exprs_data_norm)
cor_matrix <- cor(exprs_data_norm)
corrplot(cor_matrix, method = "color")
boxplot(exprs_data_norm, outline = FALSE)

dim(exprs_data_norm)
pData(gse) 
fData(gse) 
summary(exprs_data_norm)

exprs_data_log <- log2(exprs_data_norm + 1) # Log2 transformation
exprs_df_log <- as.data.frame(exprs_data_log)
cor_matrix_log <- cor(exprs_data_log)
corrplot(cor_matrix_log, method = "color")
boxplot(exprs_data_log, outline = FALSE)

dim(exprs_data_log)
summary(exprs_data_log)

group1 <- exprs_data_log[, 1:56]
group2 <- exprs_data_log[, 57:ncol(exprs_data_log)]
group3 <-exprs_data_log[, 1:25]
group4 <-exprs_data_log[, 26:50]
group5 <-exprs_data_log[, 51:76]
group6 <-exprs_data_log[, 76:ncol(exprs_data_log)]

result1 <- t.test(group1, group2)
result2 <- t.test(group3, group4)
result3 <- t.test(group5,group6)

(result1)
(result2)
(result3)

log_fc <- colMeans(group2) - colMeans(group1)
log_fc_df <- data.frame(log_fc)
log_fc

pvalues <- c(result1$p.value,result1$p.value,result1$p.value)
holm_pvalues <- p.adjust(pvalues, method = "holm")
holm_pvalues

ggplot(log_fc_df, aes(x = log_fc)) +
  geom_histogram(bins = 50) +
  labs(title = "Log Fold Change Distribution", x = "Log Fold Change")

volcano_data <- data.frame(log_fc, pvalues[1])
ggplot(volcano_data, aes(x = log_fc, y = -log10(pvalues[1]))) +
  geom_point(aes(colour = pvalues[1] < 0.05), size = 1.5) +
  scale_colour_manual(values = c("red", "black")) +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10(adjusted p-value)")

library(limma)
design <- model.matrix(~ 0 + factor(c(rep(1, 84), rep(2, ncol(exprs_data_log) - 84))))
colnames(design) <- c("group1", "group2")
fit <- lmFit(exprs_data_log, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = 2, number = nrow(exprs_data_log), adjust = "fdr")

log_fc <- results$logFC
p_values <- results$P.Value
holm_p_values <- p.adjust(p_values, method = "holm")
volcano_df <- data.frame(log_fc, -log10(holm_p_values))
volcano_df$significant <- results$adj.P.Val < 0.05

ggplot(volcano_df, aes(x = log_fc, y = -log10(holm_p_values))) +
  geom_point() +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10(p-value)") +
  scale_color_manual(values = c("red", "black"))

library(fgsea)
log_fc <- results$logFC
names(log_fc) <- rownames(results)
gene_list <- na.omit(log_fc)
gene_list <- sort(gene_list, decreasing = TRUE)
data(geneList, package = "DOSE")
gene_sets <- lapply(geneList, as.character)
gene_list <- as.numeric(gene_list)
names(gene_list) <- rownames(results)
gsea_res <- fgsea(pathways = gene_sets, stats = gene_list, nperm = 1000)
gsea_res

ego <- enrichGO(gsea_res, OrgDb="org.Hs.eg.db", keyType="ENSEMBL", ont="BP", pAdjustMethod="BH")
plotEnrichment(ego, top_term=20, showCategory=10, measure="log10p", title="GO Biological Processes")

