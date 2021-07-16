#Getting data from TCGA

#Loading essential packages
pacman::p_load(pacman, TCGAbiolinks, SummarizedExperiment, DESeq2, IHW, apeglm, pheatmap, RColorBrewer, PCAtools, reshape2)
#TCGAbiolinks,Summarized Experiment, DESeq2, IHW, biomaRT, apeglm, PCAtools - bioconductor packages
# pheatmap, RColorBrewer, reshape2 - CRAN package

#Installing essential packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#For packages that couldn't be installed using pacman
BiocManager::install("TCGAbiolinks", force = TRUE) 
BiocManager::install("biomaRt", force = TRUE)
BiocManager::install("apeglm", force = TRUE)
BiocManager::install('grimbough/biomaRt')
BiocManager::install("airway", force = T)
BiocManager::install("biobroom", force = T)


install.packages("devtools")
# Using devtools to install the package annotables
devtools::install_github("stephenturner/annotables")

options("install.lock"=FALSE)
install.packages("RcppEigen")

library(biomaRt)
 #To see loaded packages
(.packages())

#Downloading Data from TCGA 
?TCGAbiolinks
query_TCGA <-  GDCquery(
  project = "TCGA-BLCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts"
)

GDCdownload(query = query_TCGA)

getwd()
setwd("~/Desktop/Intern/")

datab <- GDCprepare(query = query_TCGA, save = T, save.filename = "data.rda")

?SummarizedExperiment
#RNA Data
rna <- as.data.frame(SummarizedExperiment::assay(datab))
save(rna, file="rna.rda")
#Clinical Data
clinical <- data.frame(datab@colData)
save(clinical, file = "clinical.rda")

table(clinical$definition)
library(psych)
describe(rna)

colnames(clinical)

#The gene names are the same in both the datasets
all(colnames(rna) == rownames(clinical))

clinical$definition
#It is recommended (but not required) to use 
#only letters, numbers, and delimiters '_' or '.', as these are safe characters for column names in R.
clinical$definition <- gsub(" ", "_", clinical$definition)

#Setting definition as factor as it will be the classification factor
clinical$definition <- as.factor(clinical$definition)

#Setting the levels
clinical$definition <- relevel(clinical$definition, ref = "Solid_Tissue_Normal")

#DESeq2
dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = clinical,
                              design = ~definition)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
View(normalized_counts)

vsd <- vst(dds, blind = TRUE)

#Extracting vst matrix
vsd_mat <- assay(vsd)

#Compute pairwise correlation analysis
vsd_cor <- cor(vsd_mat)
View(vsd_cor)

#Hierarchical clustering with correlation heatmaps

library(dplyr)
?pheatmap
pheatmap(vsd_cor, annotation = select(clinical, definition))

#Principal Component Analysis
plotPCA(vsd, intgroup="definition")
#PC1= 21%; PC2= 12%
#ToDo: Understand the reason for the small variation
#Trying to see effect of smoking on cancer
plotPCA(vsd, intgroup="paper_Number.pack.years.smoked")

plotPCA(vsd, intgroup="paper_Tobacco.smoking.history")
#Interpretation- Reformed smoker still at risk of bladder cancer?

#PC object
p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)
pairsplot(p,
          components = getComponents(p, c(1:10)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.4,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'definition',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

#Hierarchical clustering
normal_idx <- substr(colnames(assay(vsd)),14,14) == "1"
n_sample <- assay(vsd)[, c(normal_idx) ]
colnames(n_sample) <- paste("NT_", substr(colnames(n_sample),1,12))

# Dissimilarity matrix calculation
sampleDists <- dist(t(n_sample))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- c(colnames(n_sample), colnames(t_sample))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap visualization for some samples
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#############
#?plot
#plot(clinical$paper_Tobacco.smoking.history)     #Unable to plot chart
#test1 <- na.omit(clinical$paper_Tobacco.smoking.history)
#is.na(test1)
#plot(test1)
#############

#Running analysis
dds <- DESeq(dds)

# Plotting dispersion estimates
plotDispEsts(dds)

#Results
#results(dds, alpha = 0.05, altHypothesis = "greaterAbs", lfcThreshold = 1.5)
#lfcthreshold is the log2 foldchange threshold, usually useful when dealing with large number of DE genes being expressed. 

#Adding contrasts
dds_res <- results(dds,
                   contrast = c("definition", "Primary_solid_Tumor", "Solid_Tissue_Normal"),
                   alpha = 0.05,
                   altHypothesis = "greaterAbs", lfcThreshold = 1.5)
dds_res
plotMA(dds_res, ylim=c(-8,8))

#LFC Shrinkage
dds_res <- lfcShrink(dds,
                     coef=resultsNames(dds)[2], 
                     type="apeglm")
dds_res.Ordered <-  dds_res[with(dds_res, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]

plotMA(dds_res, ylim=c(-8,8)) #Shrinkage should allow for more accurate fold changes

mcols(dds_res)
head(dds_res, n=10)
summary(dds_res)


library(annotables)


# converting Ensebl id to Gene symbols using biomart
#ens2symbol<-function(ids){
#  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#  genes <- getBM(filters= "ensembl_gene_id", 
#                 attributes= c("ensembl_gene_id","hgnc_symbol"),
#                 values=ids, mart= mart)
#  return(genes)
#}
#df <- ens2symbol(row.names(dds_res))


#Human build 38 (grch38)
grch38

View(dds_res)
head(dds_res)

library(tibble)
df <- data.frame(dds_res) %>% rownames_to_column(var = "ensgene")
View(df)
head(df)


dds_res_all <- df %>%
  left_join(x = df,
            y = grch38[, c("ensgene",
                           "symbol",
                           "description")],
            by = "ensgene")
View(dds_res_all)

#Significant DEgenes -arrange
dds_res_sig <- subset(dds_res_all, padj < 0.05)
dds_res_sig <- dds_res_sig %>% arrange(padj)
View(dds_res_all)
View(dds_res_sig)

#Saving the results
write.csv(dds_res_all, file=paste0(resultsNames(dds)[2], ".csv")) #All
write.csv(dds_res_sig, file="significant.csv")

#result with Independent hypothesis weighting
resIHW <- results(dds, filterFun=ihw, alpha = 0.05, altHypothesis = "greaterAbs", lfcThreshold = 1.5)
resIHW_df <- as.data.frame(resIHW)                 
resIHW_df$ensembl_gene_id <- row.names(resIHW_df)
resIHW_df <- merge(df,resIHW_df, by = "ensembl_gene_id")

resIHWOrdered <- resIHW_df[with(resIHW_df, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]

write.csv(resIHW_df, 
          file= paste0("IHW",resultsNames(dds)[2], ".csv"))

#Visualizing results - Expression heatmap
#Subset normalized counts to significant genes

sign_norm_counts <- normalized_counts[dds_res_sig$ensgene,]


heat_colors <- brewer.pal(6, "YlOrRd")
display.brewer.all()
pheatmap(sign_norm_counts, color = heat_colors, cluster_rows = T, show_rownames = F, annotation = select(clinical, definition), scale = "row")


#Obtain logical vector regarding whether padj values are less than 0.05
dds_res_all <- dds_res_all %>% mutate(threshold = padj < 0.05)


#Volcano plot
ggplot(dds_res_all) + geom_point(aes(x = log2FoldChange, y = -log10(padj),
                                 color = threshold)) + 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

View(sign_norm_counts)
top_20 <- data.frame(sign_norm_counts)[1:20, ] %>% rownames_to_column(var="ensgene")
?gather
library(tidyr)
#top_20 <- gather(top_20, key="definition", value="normalized_counts", 2:8)
#View(top_20)
#View(clinical)
#top_20 <- inner_join(top_20, rownames_to_column(clinical, var="definition"), by="definition")
#ggplot(top_20) + geom_point(aes(x=ensgene, y=normalized_counts, color=definition)) + scale_y_log10() + xlab("Genes") + ylab("Normalized Counts") + ggtitle("Top 20 Significant DEgenes") + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust = 1)) + theme(plot.title = element_text(hjust=0.5))