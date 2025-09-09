# Analysis for gut microbial β-glucuronidases in Colorectal Cancer
Development of colorectal cancer (CRC) is accompanied by microbial and metabolic dysbiosis, with gut microbial β-glucuronidases (gmGUSs) potentially impacting carcinogenesis through de-glucuronidation of diverse important molecules. This study constructed an atlas of 550 gmGUSs from a public CRC cohort, employing 114 reference GUSs, three GUS domains, and seven conserved residues. Analysis along CRC stages revealed enrichment of Mini-Loop2 and GUS-harboring species Bacteroides cellulosilyticus and Bacteroides nordii, in late-stage. Moreover, 38 differential gmGUSs were identified totally, effectively classifying patients from controls (AUCs > 0.8). A GUSscore model based on five gmGUSs mainly from B. cellulosilyticus well predicted CRC outcomes (AUCs > 0.8). Notable gmGUS-associated genus-level shifts included reduced Lachnospira and Bifidobacterium, and increased Prevotella, Alistipes, and Fusobacterium. Particularly, several functional species were enriched in late-stage, including sulfate reducers, mucin and flavonoid degraders. Orthology-gmGUS-metabolite interactions revealed specific biological links in amino acid metabolism, vitamin biosynthesis, bacterial behavior, and LPS biosynthesis. These findings firstly define the disturbance of microbe-gmGUS-metabolite axis in colorectal tumorigenesis and its potential as early diagnostic biomarkers and therapeutic targets for CRC.

## Repository Structure and Workflow
This directory contains R scripts for reproducing all analyses and figures from the study. The scripts are organized in numerical order to facilitate sequential execution:
* 01.rarefaction_barplot.R: </br>
Generates rarefaction curves for gmGUS diversity and bar plots for reference sequences and identified gmGUSs (Figure 1b, Supplementary Figures 1-2)
* 02.LoopAnalysis.R: </br>
Analyzes loop category distribution, abundance differences, and taxonomic composition (Figure 2a-c, Supplementary Figure 3)
* 03.DiversityAnalysis.R: </br>
Performs beta-diversity analysis (PCoA) and compares alpha diversity metrics across CRC stages (Figure 2d-e)
* 04.speciesAnalysis.R: </br>
Analyzes species-level cumulative GUS abundance/number and copy number variation (CNV) (Figure 2f-h, Supplementary Figure 4)
* 05.GUSanalysis.R: </br>
Performs differential abundance analysis of gmGUSs across CRC stages and validation in independent cohorts (Figure 3a, Supplementary Figures 5, 7b)
* 06.RFmodel.R: </br>
Implements Random Forest classifiers for CRC/adenoma classification with feature selection and validation (Figure 3b-c, Supplementary Figures 6-8)
* 07.GUSscoreModel.R: </br>
Cox regression analysis and constructs GUSscore model for predicting CRC survival outcomes using LASSO (Figure 3d-g, Supplementary Figures 9-10)
* 08.speciesCorr.R: </br>
Analyzes correlations between gmGUSs and bacterial species (Figure 4, Supplementary Figure 11)
* 09.metaboliteAnalysis.R: </br>
Analyzes correlations between gmGUSs and metabolites/KEGG Orthology terms (Figure 5)
* 10.EXPandRNAseq.R: </br>
Analyzes experimental data including enzyme assays, cell experiments, and RNA-seq data (Figure 6, Supplementary Figures 12-14)

### Contents in 00.rawdata
This folder contains the raw files needed for analysis
* [CNV](00.rawdata/CNV/): Output from MIDAS2, used for copy number variation analysis.
* [Cohorts](00.rawdata/Cohorts/): The group information and GUS abundance profile of the AUS, FRA, GER cohorts.
* [FunctionalEnrichment](00.rawdata/FunctionalEnrichment/): The functional enrichment results for metabolites.
* [KOinfo](00.rawdata/KOinfo/): The KEGG-Orthology (KO) information.
* [Liter](00.rawdata/Liter/): Manual categorization of species, KOs, and metabolites.
* [RNAseq](00.rawdata/RNAseq/): Resources from RNAseq data analyses.
* [InVitro](00.rawdata/InVitro/): Data from in vitro enzyme assays and cell experiments.
* [mOTUs4](00.rawdata/mOTUs4/): Species profile based on mOTUs4.
* [supp](00.rawdata/supp/): Supplementary files from the study cohort, including profiles of species, KOs, and metabolites, as well as the differences of these features.
* [38gmGUS.tre](00.rawdata/38gmGUS.tre): The netwick tree file of 38 significant gmGUSs.
* [GUSabun_ABS.csv](00.rawdata/GUSabun_ABS.csv): The absolute abundance profile of gmGUSs.
* [GUSabun_TPM.csv](00.rawdata/GUSabun_TPM.csv): The TPM abundance profile of gmGUSs.
* [GUSsStat.csv](00.rawdata/GUSsStat.csv): The loop category, taxonomic annotation of gmGUSs.
* [GlcA.LCA.output.m8](00.rawdata/GlcA.LCA.output.m8): List of GlcA-utilizing species, tabular format.
* [clinical_supp.xlsx](00.rawdata/clinical_supp.xlsx): Clinical information of 46 samples with survival status.
* [group.csv](00.rawdata/group.csv): Group information of the study cohort.
* [referenceStat.csv](00.rawdata/referenceStat.csv): The loop category and taxonomic information of 114 references.
* [signGUSs.RDS](00.rawdata/signGUSs.RDS): A RDS file for the 38 signficant gmGUSs.


## System Requirements
The R scripts requires only a standard computer with enough RAM to support the in-memory operations.</br>
The scripts have been tested on macOS system, but Linux is theoretically feasible as well.</br>
* [R](https://cran.r-project.org) (version 4.3.2 or higher recommended)
* Required R packages:</br>
    vegan, dplyr, ggplot2, ggpubr, viridis, reshape2, ape</br>
    VennDiagram, ggtree, treeio, pheatmap, psych, coin</br>
    caret, pROC, Boruta, randomForest</br>
    survival, survminer, glmnet, timeROC, readxl</br>
    clusterProfiler, org.Hs.eg.db, grid, ggsci, ggsignif</br>
    ggbreak, Rmisc, ggprism, EnhancedVolcano, tibble</br>
* Please see each R script for its required R packages.

## Installation
1. Install R from [the official website](https://www.r-project.org)</br>
2. Install required R packages using the following command in R:</br>
`install.packages(c("vegan", "dplyr", "ggplot2", "ggpubr", "viridis", "reshape2", 
                   "ape", "VennDiagram", "ggtree", "treeio", "pheatmap", "psych", 
                   "coin", "caret", "pROC", "Boruta", "randomForest", "survival", 
                   "survminer", "glmnet", "timeROC", "readxl", "clusterProfiler", 
                   "org.Hs.eg.db", "grid", "ggsci", "ggsignif", "ggbreak", "Rmisc", 
                   "ggprism", "EnhancedVolcano", "tibble"))`</br>
if failed, try:</br>
`BiocManager::install('packageName')`</br>
3. Download or clone this repository
4. Ensure all raw data files are placed in the 00.rawdata/ directory as specified above
5. Set your working directory to the repository location:</br>
`setwd("path/to/this/directory")`</br>
6. run each script individually in R command line

## Random Forest Classifier Implementation
The Random Forest classifier (script 06.RFmodel.R) includes:
1. Feature selection using Boruta algorithm
2. Model training with hyperparameter tuning
3. Model evaluation using ROC analysis and variable importance
4. External validation in three independent cohorts (AUS, FRA, GER)

To adapt the classifier to new datasets:
∙ Ensure your data follows the same format as GUSabun_TPM.csv(rows = features, columns = samples)
∙ Provide sample metadata in the same format as group.csv
∙ Modify the feature selection parameters if needed for your dataset

## Citation
For usage of the tool, please cite the associated manuscript.

## Support
For questions regarding the code or reproducibility, please contact Junru Chen (junru.chen2019@hotmail.com).
