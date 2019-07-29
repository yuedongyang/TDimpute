# TDimpute 
TDimpute is designed to impute missing omics data by using transfer learning based deep neural network.
The method is still on progress. For any questions about the code, please contact zhoux85@mail2.sysu.edu.cn

# Requirements
tensorflow, python 3

# Data preparation
RNA-seq data (UNC IlluminaHiSeq_RNASeqV2), DNA methylation data (JHU-USC HumanMethylation450), downloaded from TCGA.

# Usage
To run script and sample dataset:

python TDimpute_without_transfer.py GPU_index cancer_name full_dataset_path imputed_dataset_path

In the script TDimpute.py, RNA_DNA_combine.csv is a 33-cancer dataset including gene expression and DNA methylation data, and not uploaded due to size limitation.

