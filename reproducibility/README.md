This folder contains the code for the downstream analyses in the manuscript.

## Dependencies

Running these analyses requires prior installation of the R libraries:
`psych` 
`PRROC`
`survival`
`plyr`
`caret`

## Supporting data

The datasets used in the analyses can be downloaded from Synapse with ID:  [syn21438134](https://www.synapse.org/#!Synapse:syn21438134). 
The `clinical_cancers` folder includes the clinical information.
The `bootstrap_cancer_V1` folder includes the actual datasets with complete gene expression data and DNA methylation data. Bootstrapping strategy is needed to sample the datasets for splitting training and testing datasets.
The `baseline_prognostic_biomarker` folder includes the prognosis-related gene list downloaded from The Human Protein Atlas.
The `WT_data` folder includes the Wilams tumor datasets from the TARGET project. 
The `UCEC_data` folder includes the UCEC dataset from the TCGA project. The samples missing gene expression are filled with NA.

## Scripts 

Here is a brief summary of what every script does:  
./Figure2: Seven imputation methods mentioned in the manuscript for generating the imputed datasets of the 16 cancers.  
./Figure3_Table2/CpG-gene_baseline.R: compute baseline correlation for the CpG-gene pairs.  
./Figure3_Table2/CpG-gene_compare.R: compute correlation of the CpG-gene pairs from datasets imputed by different methods, Figure 3.  
./Figure3_Table2/CpG-gene_select.R: compute gene lists from datasets imputed by different methods, Table 2.  
./Table3_Table4/cox_filter_baseline_generate.R: compute baseline prognosis-related gene lists.  
./Table3_Table4/cox_filter_by_cancer.R: compare the prognosis-related gene lists from different imputed datasets, Table 3.  
./Table3_Table4/cox_biomarker_by_cancer.R: compare the prognosis-related gene lists from different imputed datasets, Table 4.  
./Figure4/cindex_by_cancer.R: cluster analyses in Fig 4A.  
./Figure4/cox_filter_kmeans_by_cancer.R: survival analyses in Fig 4B.  
./UCEC and ./WT: scripts for imputing and survival analyses for dataset UCEC and WT.

