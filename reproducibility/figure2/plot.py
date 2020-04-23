import pandas as pd
import numpy as np
import seaborn as sns
sns.set()
import matplotlib.pyplot as plt

reduced_names = ['LUSC', 'KIRC', 'CESC', 'STAD', 'SARC', 'COAD','KIRP', 'LUAD', 'BLCA', 'BRCA','HNSC','LGG_','PRAD','THCA','SKCM', 'LIHC']
summary_list = np.load('/data/RMSE_list.npy')

methods = 7
RMSE_list = np.zeros([16, methods, 5, 5])
RMSE_summary = np.zeros([16, methods, 5])
RMSE_list = summary_list

cancer_c = 0
for cancertype in reduced_names:
    print(cancertype)
    RMSE_summary[cancer_c, 0, :] = np.array([np.mean(RMSE_list[cancer_c, 0, i, :]) for i in range(0, 5)]) # TDimpute
    RMSE_summary[cancer_c, 1, :] = np.array([np.mean(RMSE_list[cancer_c, 1, i, :]) for i in range(0, 5)]) # TDimpute-self
    RMSE_summary[cancer_c, 2, :] = np.array([np.mean(RMSE_list[cancer_c, 2, i, :]) for i in range(0, 5)]) # SVD
    RMSE_summary[cancer_c, 3, :] = np.array([np.mean(RMSE_list[cancer_c, 3, i, :]) for i in range(0, 5)]) # Mean method
    RMSE_summary[cancer_c, 4, :] = np.array([np.mean(RMSE_list[cancer_c, 4, i, :]) for i in range(0, 5)]) # Lasso
    RMSE_summary[cancer_c, 5, :] = np.array([np.mean(RMSE_list[cancer_c, 5, i, :]) for i in range(0, 5)]) # TOBMI
    RMSE_summary[cancer_c, 6, :] = np.array([np.mean(RMSE_list[cancer_c, 6, i, :]) for i in range(0, 5)]) # TDimpute-noTF
    np.set_printoptions(precision=3)
    print(RMSE_summary[cancer_c, :, :])
    cancer_c = cancer_c + 1

average_all_cancer = np.zeros([methods, 5])
for method_c in range(0, methods):
    for perc_c in range(0, 5):
        average_all_cancer[method_c, perc_c] = np.mean(RMSE_summary[:, method_c, perc_c])

np.set_printoptions(precision=3)
print('average_all_cancer on all cancers')
print(average_all_cancer)


dataset_all = pd.DataFrame()
cancer_c=0
samples=5
average_by_cancer = {}
std_by_cancer = {}
index_std=0
index_c = 0

for cancertype in reduced_names:
    perc_c = 0
    averages = np.zeros([methods, 5])
    stds = np.zeros([methods, 5])
    for missing_perc in [0.1, 0.3, 0.5, 0.7, 0.9]:
        DF0 = pd.DataFrame({'cancertype':[cancertype]* samples, 'Missing rate': [str(int(missing_perc * 100)) + "%"] * samples, 'method': ['TDimpute'] * samples, 'sample_count':np.arange(1,samples+1), 'RMSE': summary_list[cancer_c, 0, perc_c, :] })
        DF1 = pd.DataFrame({'cancertype':[cancertype]* samples, 'Missing rate': [str(int(missing_perc * 100)) + "%"] * samples, 'method': ['TDimpute-self'] * samples, 'sample_count':np.arange(1,samples+1), 'RMSE': summary_list[cancer_c, 1, perc_c, :] })
        DF2 = pd.DataFrame({'cancertype':[cancertype]* samples, 'Missing rate': [str(int(missing_perc * 100)) + "%"] * samples, 'method': ['SVD'] * samples, 'sample_count':np.arange(1,samples+1), 'RMSE': summary_list[cancer_c, 2, perc_c, :] })
        DF3 = pd.DataFrame({'cancertype':[cancertype]* samples, 'Missing rate': [str(int(missing_perc * 100)) + "%"] * samples, 'method': ['Mean'] * samples, 'sample_count':np.arange(1,samples+1), 'RMSE': summary_list[cancer_c, 3, perc_c, :] })
        DF4 = pd.DataFrame(
            {'cancertype': [cancertype] * samples, 'Missing rate': [str(int(missing_perc * 100)) + "%"] * samples,
             'method': ['Lasso'] * samples, 'sample_count': np.arange(1, samples + 1),
             'RMSE': summary_list[cancer_c, 4, perc_c, :]})
        DF5 = pd.DataFrame({'cancertype':[cancertype]* samples, 'Missing rate': [str(int(missing_perc * 100)) + "%"] * samples, 'method': ['TOBMI'] * samples, 'sample_count':np.arange(1,samples+1), 'RMSE': summary_list[cancer_c, 5, perc_c, :] })
        DF6 = pd.DataFrame(
            {'cancertype': [cancertype] * samples, 'Missing rate': [str(int(missing_perc * 100)) + "%"] * samples,
             'method': ['TDimpute-noTF'] * samples, 'sample_count': np.arange(1, samples + 1),
             'RMSE': summary_list[cancer_c, 6, perc_c, :]})
        dataset_all = pd.concat([dataset_all, DF0, DF1, DF2,DF3, DF4, DF5, DF6], ignore_index=True)

        method_c = 0
        for method in ['Mean', 'SVD', 'TOBMI','Lasso', 'TDimpute-noTF', 'TDimpute-self', 'TDimpute']:
            averages[method_c, perc_c] = np.mean( dataset_all[(dataset_all['cancertype'] == cancertype) & (dataset_all['Missing rate'] == str(int(missing_perc * 100)) + "%") & (dataset_all['method'] == method)]['RMSE'] )
            stds[method_c, perc_c] = np.std( dataset_all[(dataset_all['cancertype'] == cancertype) & (dataset_all['Missing rate'] == str(int(missing_perc * 100)) + "%") & (dataset_all['method'] == method)]['RMSE'] )
            method_c = method_c + 1
        perc_c = perc_c + 1
        average_by_cancer[cancertype] = averages
        std_by_cancer[cancertype] = stds

    cancer_c = cancer_c + 1

#############
dataset_mean_cancer = pd.DataFrame()
average_all_cancer = np.zeros([7, 5])
std_all_cancer = np.zeros([7, 5])
perc_c = 0
index_c = 0
for missing_perc in ['10%','30%','50%','70%','90%']:
    method_c = 0
    for method in ['Mean','SVD', 'TOBMI','Lasso','TDimpute-noTF', 'TDimpute-self', 'TDimpute']:
        for sample_count in [1, 2, 3, 4, 5]:
            ##mean of 16 cancers, first step   !!!!
            aa = np.mean(dataset_all[(dataset_all['sample_count'] == sample_count) & (dataset_all['method'] == method) & (dataset_all['Missing rate'] == missing_perc)]['RMSE'])
            DF = pd.DataFrame({'Missing rate': missing_perc, 'method': method, 'sample_count': sample_count, 'RMSE': aa}, index=[index_c])
            dataset_mean_cancer = pd.concat([dataset_mean_cancer, DF], ignore_index=True)
            index_c=index_c+1

        # mean of 5 samples, second step !!!!
        average_all_cancer[method_c, perc_c] = np.mean( dataset_mean_cancer[(dataset_mean_cancer['method'] == method) & (dataset_mean_cancer['Missing rate'] == missing_perc)]['RMSE'] )
        std_all_cancer[method_c, perc_c] = np.std( dataset_mean_cancer[(dataset_mean_cancer['method'] == method) & (dataset_mean_cancer['Missing rate'] == missing_perc)]['RMSE'] )
        method_c = method_c + 1
    perc_c = perc_c + 1

print('average_all_cancer on all cancers')
print(average_all_cancer)

##########average by 16 cancers
fig, (ax) = plt.subplots(1,1)
x = np.array(['10%','30%','50%','70%','90%'])
y = average_all_cancer
yerr = std_all_cancer
linewidth=2
ax.errorbar(x, y[0,:], yerr=yerr[0,:], label='Mean', linestyle=(0, (6, 1, 2, 1,2,1)), linewidth=linewidth, color = 'r')
ax.errorbar(x, y[1,:], yerr=yerr[1,:], label='SVD', linestyle=":", linewidth=linewidth, color = 'g')
ax.errorbar(x, y[2,:], yerr=yerr[2,:], label='TOBMI', linestyle=(0, (3, 1, 1, 1,1,1)), linewidth=linewidth, color = 'orangered')
ax.errorbar(x, y[3,:], yerr=yerr[3,:], label='Lasso', linestyle=(0, (3, 1, 4, 1,1,1)), linewidth=linewidth, color = 'purple')
ax.errorbar(x, y[4,:], yerr=yerr[4,:], label='TDimpute-noTF', linestyle="-.", linewidth=linewidth, color = 'cornflowerblue')
ax.errorbar(x, y[5,:], yerr=yerr[5,:], label='TDimpute-self', linestyle="--", linewidth=linewidth, color = 'b')
ax.errorbar(x, y[6,:], yerr=yerr[6,:], label='TDimpute', linestyle="-", linewidth=linewidth, color = 'royalblue')

# get handles
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] for h in handles]
# use them in the legend
plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=7)
ax.set_ylabel('Mean absolute error')
ax.set_xlabel('Missing rate')