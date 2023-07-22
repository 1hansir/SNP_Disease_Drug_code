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

tissue_dic = {'GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt': 'breast',
              'GTEx/GTEx_Analysis_v8_eQTL/Colon_Sigmoid.signifpairs.txt' : 'colon',
            'GTEx/GTEx_Analysis_v8_eQTL/Kidney_Cortex.signifpairs.txt': 'kidney',
              'GTEx/GTEx_Analysis_v8_eQTL/Liver.signifpairs.txt': 'liver',
             'GTEx/GTEx_Analysis_v8_eQTL/Lung.signifpairs.txt':'lung',
             'GTEx/GTEx_Analysis_v8_eQTL/Prostate.signifpairs.txt' : 'prostate',
             'GTEx/GTEx_Analysis_v8_eQTL/Skin_Not_Sun_Exposed_Suprapubic.signifpairs.txt' : 'skin'
            }

gene_index_filename = 'mapping_files/Gene_INDEX/Gene_index.csv'
snp_index_filename = 'mapping_files/SNP_INDEX/snpid_file.csv'
gene_index = pd.read_csv(mdir + gene_index_filename, sep='\t')
snp_index = pd.read_csv(mdir + snp_index_filename, sep = '\t', low_memory=False)
gene_index_dic = dict(zip(gene_index['gene_id'],gene_index['gene_index']))
snp_index_dic = dict(zip(snp_index['variant_id'],snp_index['SNP_index']))

# read gtex data from different tissues:
gtex_alldata = []
for tissue_file in gtex_filenames:
    gtex_data_t = pd.read_csv(mdir + tissue_file, sep='\t')
    if tissue_file == 'GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt':
        gtex_alldata = gtex_data_t
    else:
        gtex_alldata = pd.concat([gtex_alldata, gtex_data_t])
print(gtex_alldata.columns)
all_pairs = gtex_alldata[['variant_id','gene_id','slope']]
all_pairs['gene_id'] = all_pairs['gene_id'].str.split('.').str[0]
all_pairs = all_pairs[all_pairs['variant_id'].isin(snp_index_dic.keys())]

all_pairs['gene_index'] = [gene_index_dic[x] for x in all_pairs['gene_id'] ]
all_pairs['snp_index'] = [snp_index_dic[x] for x in all_pairs['variant_id'] ]
all_pairs.drop_duplicates(subset = ['gene_index','snp_index'],keep = 'first', inplace = True)
all_pairs = all_pairs.reset_index(drop=True)

all_pairs.to_csv(mdir + 'mapping_files/PAIRS/all_pairs_t.csv', index_label=['pair_label'])

'''  SLOW algorithm
for i in range(0, len(all_pairs)):
    gene_id_i = all_pairs.loc[i, 'gene_id']
    gene_index_i = gene_index_dic[gene_id_i]
    all_pairs.loc[i, 'gene_index'] = int(gene_index_i)

    snp_id_i = all_pairs.loc[i, 'variant_id']
    try:
        snp_index_i = snp_index_dic[snp_id_i]
    except KeyError:
        snp_index_i = -1
    all_pairs.loc[i, 'snp_index'] = int(snp_index_i)

    if i % 10000 == 0:
        print(i / len(all_pairs))

all_pairs = all_pairs[all_pairs['snp_index'] != -1]
'''

for interest_tissue in gtex_filenames:
    tissue_pairs = pd.read_csv(mdir + interest_tissue, sep='\t',usecols=['variant_id', 'gene_id','slope'])
    tissue_pairs['gene_id'] = tissue_pairs['gene_id'].str.split('.').str[0]

    tissue_pairs = tissue_pairs[tissue_pairs['variant_id'].isin(snp_index_dic.keys())]
    tissue_pairs = tissue_pairs.reset_index(drop=True)

    tissue_pairs['gene_index'] = [gene_index_dic[x] for x in tissue_pairs['gene_id'] ]
    tissue_pairs['snp_index'] = [snp_index_dic[x] for x in tissue_pairs['variant_id']]

    tissue_pairs.to_csv(mdir + 'mapping_files/PAIRS/tissue_{}_t.csv'.format(tissue_dic[interest_tissue]), index_label=['pair_index'])

