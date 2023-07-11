import pandas as pd


mdir = 'D:/Pycharm_projects/SNP_Disease_Drug/'

cmap_filename = mdir + 'CMap/Format_data/cp_MCF7_lm.csv'
cmap_alldata = pd.read_csv(cmap_filename)
lm_gene_ls = cmap_alldata.columns[1:]

gtex_filename = mdir + 'GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt'

gtex_alldata = pd.read_csv(gtex_filename,sep='\t')

all_slope = gtex_alldata['slope']
all_SNPs = gtex_alldata['variant_id']
all_GENEs = gtex_alldata['gene_id']
n_gtex = len(all_GENEs)

uni_GENEs = all_GENEs.unique()

# create a dictionary to restore SNP<-->Gene pairs
pair_dic = {}
pair_dic_num = {}
for i in range(0,n_gtex):
    gene_id = all_GENEs[i]
    if gene_id not in pair_dic.keys():
        pair_dic[gene_id] = []
        pair_dic_num[gene_id] = 0

    pair_dic[gene_id].append(all_SNPs[i] + '-' + str(all_slope[i]))
    pair_dic_num[gene_id] += 1

# SNP-centered pair dictionary
snp_pair_dic = {}
snp_pair_dic_num = {}
for i in range(0,n_gtex):
    snp_id = all_SNPs[i]
    if snp_id not in snp_pair_dic.keys():
        snp_pair_dic[snp_id] = []
        snp_pair_dic_num[snp_id] = 0

    snp_pair_dic[snp_id].append(all_GENEs[i] + '-' + str(all_slope[i]))
    snp_pair_dic_num[snp_id] += 1

import pickle

with open(mdir + 'GTEx/Format_data/saved_gene_dict.pkl', 'wb') as f:
    pickle.dump(pair_dic, f)
with open(mdir + 'GTEx/Format_data/saved_snp_dict.pkl', 'wb') as f:
    pickle.dump(snp_pair_dic, f)

with open(mdir + 'GTEx/Format_data/saved_gene_num.pkl', 'wb') as f:
    pickle.dump(pair_dic_num, f)
with open(mdir + 'GTEx/Format_data/saved_snp_num.pkl', 'wb') as f:
    pickle.dump(snp_pair_dic_num, f)






