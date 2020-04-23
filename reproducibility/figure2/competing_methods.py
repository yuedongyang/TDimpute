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

cancer_names = ['BRCA', "SARC", "LUSC", "BLCA", "KIRC", "LGG_", "PRAD", "LUAD", "LIHC", "SKCM",
                 "HNSC", "CESC", "COAD", "KIRP", "THCA", "STAD"]

RNA_size = 19027
DNA_size = 27717

sample_size = 5
cancer_num = 16
loss_list_Mean = np.zeros([cancer_num, 5, sample_size])
loss_list_SVD = np.zeros([cancer_num, 5, sample_size])
loss_list_Lasso = np.zeros([cancer_num, 5, sample_size])
cancer_c = 0
for cancertype in cancer_names:
    perc = 0
    for missing_perc in [0.1,0.3,0.5,0.7,0.9]:
        for sample_count in range(1, sample_size + 1):
            shuffle_cancer = pd.read_csv(datadir+'/bootstrap_cancer_V1/'+cancertype+str(int(missing_perc*100))+'_'+str(sample_count)+'.csv',
                                         delimiter=',', index_col = 0, header = 0)          
            ## 16.5 is for gene expression normalization in [0, 1]
            aa = np.concatenate((shuffle_cancer.values[:,:RNA_size]/16.5, shuffle_cancer.values[:,RNA_size:]), axis=1 )
            shuffle_cancer = pd.DataFrame(aa, index=shuffle_cancer.index, columns = shuffle_cancer.columns)

            RDNA = shuffle_cancer.values
            test_data = RDNA[0:int(RDNA.shape[0]*missing_perc),:]
            train_data = RDNA[int(RDNA.shape[0]*missing_perc):,:]

            corrupted_holdout = test_data.copy()
            corrupted_holdout[:,:19027] = np.nan
            df_combine = pd.DataFrame(np.concatenate([corrupted_holdout, train_data], axis=0))
            print('name:', cancertype, ' missing rate:', missing_perc, 'train datasize:', train_data.shape, ' test datasize: ', test_data.shape)

			############## Mean method
            X_filled = SimpleFill(fill_method="mean").fit_transform(df_combine)
            RNA_txt = pd.DataFrame(X_filled[:, :RNA_size]* 16.5 , index=shuffle_cancer.index, columns = shuffle_cancer.columns[:RNA_size])
            RNA_txt.to_csv(datadir+'/filled_data/Mean_'+ cancertype + str(missing_perc * 100) + '_' + str(
                sample_count)+'.csv')

            nz = test_data[:,:RNA_size].size
            nnm_mse = np.sqrt((np.linalg.norm((X_filled[:test_data.shape[0],:RNA_size] - test_data[:,:RNA_size])* 16.5) ** 2) / nz)
            print("Mean method, RMSE: %f" % nnm_mse)
            loss_list_Mean[cancer_c, perc, sample_count - 1] = nnm_mse
			
			
			##############SVD
			rank = 10
            X_filled = IterativeSVD(rank, init_fill_method="mean", verbose=False,convergence_threshold=0.0000001).fit_transform(df_combine)
            RNA_txt = pd.DataFrame(X_filled[:, :RNA_size]* 16.5 , index=shuffle_cancer.index, columns = shuffle_cancer.columns[:RNA_size])
            RNA_txt.to_csv(datadir+'/filled_data/SVD_'+ cancertype + str(missing_perc * 100) + '_' + str(
                sample_count)+'.csv')

            nz = test_data[:,:RNA_size].size
            nnm_mse = np.sqrt((np.linalg.norm((X_filled[:test_data.shape[0],:RNA_size] - test_data[:,:RNA_size])* 16.5) ** 2) / nz)
            print("SVD, RMSE: %f" % nnm_mse)
            loss_list_SVD[cancer_c, perc, sample_count - 1] = nnm_mse

            ##############Lasso 
            y = train_data[:, :RNA_size]  ## gene expression
            X = train_data[:, RNA_size:]  ## DNA methylation

            starttime = datetime.datetime.now()
            reg = Lasso(alpha=0.1, random_state=0).fit(X, y)
            reconstruct = reg.predict(test_data[:, RNA_size:])
            nnm_mse = np.sqrt(mean_squared_error(test_data[:, :RNA_size]* 16.5, reconstruct* 16.5))
			print("Lasso, RMSE: %f" % nnm_mse)
            endtime = datetime.datetime.now()
            print('Time elapsed: ',(endtime - starttime).seconds)
            loss_list_Lasso[cancer_c, perc, sample_count - 1] = nnm_mse

            filled_data = np.concatenate([reconstruct, train_data[:, :RNA_size]], axis=0)
            RNA_txt = pd.DataFrame(filled_data[:, :RNA_size]* 16.5, index=shuffle_cancer.index, columns = shuffle_cancer.columns[:RNA_size])
            RNA_txt.to_csv(datadir+'/filled_data/Lasso_'+ cancertype + str(missing_perc * 100) + '_' + str(
                sample_count)+'.csv')

        perc = perc + 1
    print([np.mean(loss_list_Mean[cancer_c, i, :]) for i in range(0, 5)])
    print([np.mean(loss_list_SVD[cancer_c, i, :]) for i in range(0, 5)])
    print([np.mean(loss_list_Lasso[cancer_c, i, :]) for i in range(0, 5)])
    cancer_c = cancer_c + 1
