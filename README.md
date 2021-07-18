# Differential-Expression-Analysis
### Differential gene expression analysis of RNA-seq data from TCGA

**Code: [main.R](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/main.R)**

**Result Files:**

**-[Result](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/definition_Primary_solid_Tumor_vs_Solid_Tissue_Normal.csv)**

**-[Result with Significant Genes](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/significant.csv)**

**-[Result with Independent Hypothesis Weighting](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/IHWdefinition_Primary_solid_Tumor_vs_Solid_Tissue_Normal.csv)**

**Pairs plot**
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/Pairs%20plot.png?raw=true)

**Principal Component Analysis**
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/Principal%20Component%20Analysis(definition).png?raw=true)
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/PCA(Number%20of%20year%20smoked).png?raw=true)
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/PCA(Smoking%20history).png?raw=true)

**Heatmap visualization of some samples to get preliminary idea about the data**
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/Heatmap%20visualization-some%20samples.png?raw=true)

**Hierarchial clustering Heatmap**

Hierarchical clustering with heatmaps is used to assess the similarity in gene expression between the different samples in a dataset. This technique is used to explore how similar replicates are to each other and whether the samples belonging to different sample groups cluster separately. The heatmap is created by using the gene expression correlation values for all pairwise combinations of samples in the dataset, with the value 1 being perfect correlation. The hierarchical tree shows which samples are more similar to each other and the colors in the heatmap depict the correlation values
![Hierarchial clustering Heatmap](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/Heatmap%20clustering%20Heamap_edited.png?raw=true)


**Dispersion estimate**
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/Dispersion%20estimates.png?raw=true)

**Adding contrasts**
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/LFC%20contrasts.png?raw=true)

**LFC shrinkage**
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/LFC%20Shrinkage.png?raw=true)

**Expression heatmap - Subset normalized counts for significant genes**
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/Heatmap%20normalized%20counts%20significant%20genes.png?raw=true)

**Volcano plot**
![](https://github.com/rushil1904/Differential-Expression-Analysis/blob/main/Volcano%20plot.png?raw=true)
