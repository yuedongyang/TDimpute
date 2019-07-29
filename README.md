# TDimpute 
TDimpute is designed to impute missing gene expression data from DNA methylation data by using transfer learning based deep neural network.
The method is still on progress. For any questions about the code, please contact zhoux85@mail2.sysu.edu.cn

# Requirements
tensorflow, python 3

# Data preparation
RNA-seq data (UNC IlluminaHiSeq_RNASeqV2), DNA methylation data (JHU-USC HumanMethylation450), downloaded from TCGA.
To simulate omics missing dataset, we randomly select an increasing fraction (10%, 30%, 50%, 70%, 90%) of samples in the full dataset (gold standard) and remove their gene expression data. The samples with missing gene expression are set as testing dataset and the remaining samples with complete omics are set as training dataset. At each level of missing rate, we repeat this procedure 5 times to obtain a robust evaluation.

# Usage
To run script and sample dataset:
python TDimpute_without_transfer.py GPU_index cancer_name full_dataset_path imputed_dataset_path

In the script TDimpute.py, RNA_DNA_combine.csv is a 33-cancer dataset including gene expression and DNA methylation data, and not uploaded due to size limitation.

