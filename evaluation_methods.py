import pandas as pd
import numpy as np

np.set_printoptions(precision=3)
cancertype = 'WT'

original_dat_path_RNA = '/data0/zhoux/UCSC_RNA_WT.csv'
RNA_target = pd.read_csv(original_dat_path_RNA, delimiter=',', index_col=0, header=0)
RNA_target.index = [x[:19] for x in RNA_target.index.values]
RNA_target = RNA_target.values

missing_perc = 0.5
missing_size = int(RNA_target.shape[0] * missing_perc)
test_data = RNA_target[0:missing_size, :]
nz = test_data.shape[0] * test_data.shape[1]

imputed_RNA = pd.read_csv('TDimpute_'+cancertype+str(int(missing_perc * 100)) + '.csv', index_col=0, header=0)
imputed_RNA = imputed_RNA.iloc[0:missing_size, :].values

## Four evaluation metrics for data imputation
## RMSE (Root Mean Square Error), MAE (Mean Absolute Error), PCC (Person Correlation Coefficient), MAPE (Mean Average Percentage Error).
RMSE_test = np.sqrt(np.sum((imputed_RNA - test_data) ** 2) / nz)
MAE = np.sum(abs(imputed_RNA - test_data)) / nz
aa = abs((imputed_RNA - test_data) / test_data)
MAPE = np.mean(aa[~np.isinf(aa) & ~np.isnan(aa)])

PCC_arr = [] # PCC (by sample) between the imputed data and the original full data
for sample_pcc in range(0, missing_size):
    PCC = pearsonr(imputed_RNA[sample_pcc, :], test_data[sample_pcc, :])[0]
    PCC_arr.append(PCC ** 2)

PCC_by_gene = [] # PCC (by gene) between the imputed data and the original full data
for sample_gene in range(0, imputed_RNA.shape[1]):
    PCC = pearsonr(imputed_RNA[:, sample_gene], test_data[:missing_size, sample_gene])[0]
    PCC_by_gene.append(PCC ** 2)

print('TDimpute evaluation: ', RMSE_test, MAE, MAPE, np.mean(PCC_arr), np.nanmean(PCC_by_gene))
