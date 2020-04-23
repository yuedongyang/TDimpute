np.set_printoptions(precision=2)
import random
import pandas as pd
import numpy as np
import os
import time
import datetime
from fancyimpute import SimpleFill, IterativeSVD, Lasso
from sklearn.metrics import mean_squared_error

datadir = '/data'
DNA_target = pd.read_csv(datadir+'/quantiles_DNA_WT_RSEM.csv', delimiter=',',index_col=0, header=0)
RNA_target = pd.read_csv(datadir+'/quantiles_RNA_WT_RSEM.csv', delimiter=',', index_col=0, header=0)
DNA_target.index = [x[:19] for x in DNA_target.index.values]
shuffle_cancer = pd.merge(RNA_target, DNA_target, left_index=True, right_index=True, how = 'inner')
RNA_size = RNA_target.shape[1]
DNA_size = 23003
cancertype = 'WT'

sample_size = 5
cancer_num = 1
loss_list_Mean = np.zeros([cancer_num, 5, sample_size])
loss_list_SVD = np.zeros([cancer_num, 5, sample_size])
loss_list_Lasso = np.zeros([cancer_num, 5, sample_size])
cancer_c = 0

perc = 0
for missing_perc in [0.5]:
    for sample_count in range(1, sample_size + 1):
        ## train/test data split
        train_data = shuffle_cancer.sample(frac=(1 - missing_perc), random_state=sample_count, axis=0, replace=False)
        test_data = shuffle_cancer[~shuffle_cancer.index.isin(train_data.index)]  # bool index, e.g. df[df.A>0]
        new_dataset = pd.concat([test_data, train_data], axis=0)
        train_data = train_data.values
        test_data = test_data.values
        print('train datasize:', train_data.shape, ' test datasize: ', test_data.shape)

        corrupted_holdout = test_data.copy()
        corrupted_holdout[:,:RNA_size] = np.nan
        df_combine = pd.DataFrame(np.concatenate([corrupted_holdout, train_data], axis=0))
        print('name:', cancertype, ' missing rate:', missing_perc, 'train datasize:', train_data.shape, ' test datasize: ', test_data.shape)

        ############## Mean method
        X_filled = SimpleFill(fill_method="mean").fit_transform(df_combine)
        RNA_txt = pd.DataFrame(X_filled[:, :RNA_size] , index=shuffle_cancer.index, columns = shuffle_cancer.columns[:RNA_size])
        RNA_txt.to_csv(datadir+'/filled_data/Mean_'+ cancertype + str(missing_perc * 100) + '_' + str(
            sample_count)+'.csv')

        nz = test_data[:,:RNA_size].size
        nnm_mse = np.sqrt((np.linalg.norm((X_filled[:test_data.shape[0],:RNA_size] - test_data[:,:RNA_size])) ** 2) / nz)
        print("Mean method, RMSE: %f" % nnm_mse)
        loss_list_Mean[cancer_c, perc, sample_count - 1] = nnm_mse

        ##############SVD
        rank = 10
        X_filled = IterativeSVD(rank, init_fill_method="mean", verbose=False,convergence_threshold=0.0000001).fit_transform(df_combine)
        RNA_txt = pd.DataFrame(X_filled[:, :RNA_size] , index=shuffle_cancer.index, columns = shuffle_cancer.columns[:RNA_size])
        RNA_txt.to_csv(datadir+'/filled_data/SVD_'+ cancertype + str(missing_perc * 100) + '_' + str(
            sample_count)+'.csv')

        nz = test_data[:,:RNA_size].size
        nnm_mse = np.sqrt((np.linalg.norm((X_filled[:test_data.shape[0],:RNA_size] - test_data[:,:RNA_size])) ** 2) / nz)
        print("SVD, RMSE: %f" % nnm_mse)
        loss_list_SVD[cancer_c, perc, sample_count - 1] = nnm_mse

        ##############Lasso
        y = train_data[:, :RNA_size]  ## gene expression
        X = train_data[:, RNA_size:]  ## DNA methylation

        starttime = datetime.datetime.now()
        reg = Lasso(alpha=0.1, random_state=0).fit(X, y)
        reconstruct = reg.predict(test_data[:, RNA_size:])
        nnm_mse = np.sqrt(mean_squared_error(test_data[:, :RNA_size], reconstruct))
        print("Lasso, RMSE: %f" % nnm_mse)
        endtime = datetime.datetime.now()
        print('Time elapsed: ',(endtime - starttime).seconds)
        loss_list_Lasso[cancer_c, perc, sample_count - 1] = nnm_mse

        filled_data = np.concatenate([reconstruct, train_data[:, :RNA_size]], axis=0)
        RNA_txt = pd.DataFrame(filled_data[:, :RNA_size], index=shuffle_cancer.index, columns = shuffle_cancer.columns[:RNA_size])
        RNA_txt.to_csv(datadir+'/filled_data/Lasso_'+ cancertype + str(missing_perc * 100) + '_' + str(
            sample_count)+'.csv')

    perc = perc + 1
print([np.mean(loss_list_Mean[cancer_c, i, :]) for i in range(0, 5)])
print([np.mean(loss_list_SVD[cancer_c, i, :]) for i in range(0, 5)])
print([np.mean(loss_list_Lasso[cancer_c, i, :]) for i in range(0, 5)])
cancer_c = cancer_c + 1
