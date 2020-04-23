# TDimpute 
TDimpute is designed to impute missing gene expression data from DNA methylation data by using transfer learning based neural network. 
The preprint paper could be found at [BioRxiv](https://doi.org/10.1101/803692). For any questions about the code or original datasets, please contact zhoux85@mail2.sysu.edu.cn

# Requirements
tensorflow 1.11.0  
python 3.6.5  
preprocessCore 1.48.0

# Data preparation
 All the datasets below can be downloaded from Synapse with ID: [syn21438134](https://www.synapse.org/#!Synapse:syn21438134).

* In the experiments, we use RNA-seq data (UNC IlluminaHiSeq_RNASeqV2_RSEM) and DNA methylation data (JHU-USC HumanMethylation450) from TCGA. To reproduce the main simulation results in our paper, please use this dataset: ```RNA_DNA_combine.csv``` 

* ```quantiles_DNA_TCGA_RSEM.csv``` and ```quantiles_RNA_TCGA_RSEM.csv``` are two quantile normalized datasets, which can be used to pretrain the pancancer model on TCGA. For dataset outside TCGA, we recommend using these two percentile version data if the non-TCGA data is normalized by percentile ranking.

* ```reference_distribution_DNA_TCGA.RData``` (DNA methylation) and ```reference_distribution_RNA_TCGA.RData``` (RNA-seq) are two processed files using funciton "normalize.quantiles.determine.target" in R package "preprocessCore". They can be loaded directly as reference distribution of DNA methylation and RNA-seq data from TCGA.

* We use the Wilms tumor (WT) dataset (```DNA_WT.csv``` and ```UCSC_RNA_WT.csv```) from TARGET cancer project as a example for imputing RNA-seq data using DNA methylation data. Note that the RNA-seq data should be quantified as RSEM estimated read counts, because we pretrained the neural network with the RNA-seq data of RNASeqV2_RSEM. The pretrained model with other quantification, such as readcounts, TPM, will be provided later.

```quantile_normalization_process.R``` is used to remove technical variabilities between TCGA and the dataset (outside TCGA) you want to impute: specifically, the TCGA data is considered as reference to normalize the your dataset into the same distribution. 

# Pretrained model
```ref_general_model_quantiles.ckpt``` is the pretrained model on quantile normalized TCGA pancancer data (i.e., ```quantiles_DNA_TCGA_RSEM.csv``` and ```quantiles_RNA_TCGA_RSEM.csv```). It can be downloaded from Synapse with ID: [syn21438134](https://www.synapse.org/#!Synapse:syn21438134).

# Files
Code to reproduce the main results in the paper:

* ```TDimpute_self.py``` implements the training of the neural network from scratch to impute RNA-seq data and the evaluation of imputation performance.

* ```TDimpute.py``` implements the training of the neural network whose weights and biases are initialized based on the values learnt previously on the pancancer dataset.

* To run script on the example dataset (Wilms tumor from TARGET project):
```python TDimpute_example.py quantiles_DNA_TCGA_RSEM.csv quantiles_RNA_TCGA_RSEM.csv DNA_WT.csv UCSC_RNA_WT.csv ./```

* ```evaluation_methods.py``` Different methods to evaluate the imputation accuracy.

Example code for RNA-seq imputation:
* ```TDimpute_finetune.py``` fine-tune the pretrained model if both methylation and RNA-seq available.
* ```TDimpute_pretrain.py``` uses the pretrained model directly for imputation if no RNA-seq data are provided for fine-tune (only methylation).

# Reproducibility

Please refer to the corresponding `reproducibility` folder for detailed information.


# Citation
If you find this work useful in your research, please consider citing our paper: "Imputing missing RNA-seq data from DNA methylation by using transfer learning based neural network", [BioRxiv](https://doi.org/10.1101/803692), 2019.
