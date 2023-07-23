import pandas as pd
import numpy as np

mdir = 'D:/Pycharm_projects/SNP_Disease_Drug/'

snp_pair_filename = 'mapping_files/PAIRS/all_pairs_t.csv'
cmap_oe_filename = 'CMap/Format_data/Maxpool/oe.csv'
cmap_ko_filename = 'CMap/Format_data/Maxpool/ko.csv'
gene_cmap_info_filename = "/CMap/Datasets/geneinfo_beta.txt"
snp_pairs = pd.read_csv(mdir + snp_pair_filename)
oe_vec = pd.read_csv(mdir + cmap_oe_filename,index_col=['gene'])
oe_vec_init = pd.read_csv(mdir + cmap_oe_filename)
ko_vec = pd.read_csv(mdir + cmap_ko_filename,index_col=['gene'])
cmap_gene_file = pd.read_csv(mdir + gene_cmap_info_filename, sep="\t", dtype=str)
gene_info = cmap_gene_file['ensembl_id'].unique()


vec_0 = 0 * oe_vec_init.drop(['gene'],axis=1).iloc[0]
vec_0 = 0 * oe_vec.iloc[0]

print(vec_0)
print(snp_pairs.columns)
print(oe_vec.columns)

gene_oe = oe_vec.index.unique()
gene_ko = ko_vec.index.unique()
gene_snp = snp_pairs['gene_id'].unique()
print(len(gene_snp),len(gene_oe),len(gene_ko), len(gene_info))

gene_oe_set = set(gene_oe)
gene_ko_set = set(gene_ko)
gene_snp_set = set(gene_snp)
gene_cmap_set = set(gene_info)

print(len(gene_snp_set & gene_ko_set))
print(len(gene_snp_set & gene_oe_set))
print(len(gene_snp_set & gene_oe_set & gene_ko_set))
print(len(gene_cmap_set & gene_snp_set))

snp_list = snp_pairs['snp_index'].unique()
snp_vec_all = vec_0
flag = 0
count_miss = 0
snp_name = []

for i in range(0, round(len(snp_list)/1000)):
    snp_index_i = snp_list[i]
    pairs = snp_pairs[snp_pairs['snp_index'] == snp_index_i]
    snp_vec_i = vec_0
    vec_add = vec_0

    gene_list_i = list(pairs['gene_id'])
    slope_list_i = list(pairs['slope'])
    flag_i = 0
    for j in range(0,len(gene_list_i)):
        slope_i_j = slope_list_i[j]
        gene_i_j = gene_list_i[j]

        if slope_i_j > 0 and (gene_i_j in gene_oe):
            vec_add = slope_i_j * oe_vec.loc[gene_i_j]
            flag_i = 1
        elif slope_i_j < 0 and (gene_i_j in gene_ko):
            vec_add = slope_i_j * ko_vec.loc[gene_i_j]
            flag_i = 1
        else:
            count_miss += 1

        snp_vec_i = snp_vec_i + vec_add

    if i % 10000 == 0:
            print(i/len(snp_list))

    if flag_i == 0:
        continue

    snp_vec_all = pd.concat([snp_vec_all,snp_vec_i],axis=1)
    snp_name.append(snp_index_i)


snp_vec_all = pd.DataFrame(snp_vec_all)
snp_name.insert(0, 'index')
snp_vec_all.columns = snp_name
snp_vec_all = snp_vec_all.T.drop('index')
# snp_vec_all.insert(0,'snp_index',snp_pairs['snp_index'].unique())
snp_vec_all.to_csv(mdir + 'GTEx/Format_data/snp_vec.csv')

print('missing rate:{}'.format(count_miss/len(snp_pairs)))