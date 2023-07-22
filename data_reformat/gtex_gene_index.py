import pandas as pd
import numpy as np

mdir = 'D:/Pycharm_projects/SNP_Disease_Drug/'
gtex_filenames = ['GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Colon_Sigmoid.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Kidney_Cortex.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Liver.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Lung.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Prostate.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Skin_Not_Sun_Exposed_Suprapubic.signifpairs.txt']

# read gtex data from different tissues:
gtex_alldata = []
for tissue_file in gtex_filenames:
    gtex_data_t = pd.read_csv(mdir + tissue_file, sep='\t')
    if tissue_file == 'GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt':
        gtex_alldata = gtex_data_t
    else:
        gtex_alldata = pd.concat([gtex_alldata, gtex_data_t])

print(gtex_alldata.columns)
all_GENEs = gtex_alldata['gene_id'].unique()
all_GENEs = pd.Series(all_GENEs)
all_GENEs = pd.DataFrame(all_GENEs.str.split('.').str[0])

n_gtex = len(all_GENEs)
print("# of GENEs(unique) in GTEx file:{}".format(n_gtex))

print(all_GENEs)
all_GENEs.to_csv(mdir+'mapping_files/Gene_INDEX/Gene_index.csv',sep='\t',header=['gene_id'],index_label=['gene_index'])

