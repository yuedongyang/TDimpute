# TDimpute 
TDimpute is designed to impute missing gene expression data from DNA methylation data by using transfer learning based neural network.
The preprint paper could be found at [BioRxiv](https://doi.org/10.1101/803692). For any questions about the code or original datasets, please contact zhoux85@mail2.sysu.edu.cn

# Requirements
tensorflow 1.11.0  
python 3.6.5  
preprocessCore 1.48.0

# Data preparation
RNA-seq data (UNC IlluminaHiSeq_RNASeqV2_RSEM), DNA methylation data (JHU-USC HumanMethylation450), downloaded from TCGA.

We use the Wilms tumor dataset from TARGET cancer project as a example for imputing RNA-seq data using DNA methylation data. Note that the RNA-seq data should be quantified as RSEM estimated read counts, since we pretrained the neural network with the RNA-seq data of RNASeqV2_RSEM. The pretrained model with other quantification, such as readcounts, TPM, will be provided later.

The *pretrained model, example datasets (Wilms tumor), quantiles_DNA(RNA)_TCGA_RSEM.csv reference_distribution_DNA(RNA)_TCGA.RData* can be downloaded from Synapse with ID: [syn21438134](https://www.synapse.org/#!Synapse:syn21438134).

# Usage
### quantile normalization
quantile_normalization_process.R is used to remove technical variabilities between TCGA and the dataset (outside TCGA) you want to impute: specifically, the TCGA data is considered as reference to normalize the your dataset into the same distribution.  
"quantiles_DNA_TCGA_RSEM.csv" and "quantiles_RNA_TCGA_RSEM.csv" are two quantile normalized datasets, which can be used to pretrain the pancancer model on TCGA. For dataset outside TCGA, we recommend using these two percentile version data if the non-TCGA data is normalized by percentile ranking.
"reference_distribution_DNA_TCGA.RData" and "reference_distribution_RNA_TCGA.RData" are two processed files using funciton "normalize.quantiles.determine.target" in R package "preprocessCore". They can be loaded directly as reference distribution of DNA methylation and RNA-seq data from TCGA.

### To run script and sample dataset:
python TDimpute.py GPU_index cancer_name full_dataset_path imputed_dataset_path
