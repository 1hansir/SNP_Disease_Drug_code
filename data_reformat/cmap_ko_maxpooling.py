from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import numpy as np


gene_info = pd.read_csv("../../CMap/Datasets/geneinfo_beta.txt", sep="\t", dtype=str)
HUGO_emb_dict = dict(zip(gene_info['gene_symbol'],gene_info['ensembl_id']))
hugo_ens_sup = pd.read_csv('../../mapping_files/HUGOSymbols_enid/all_hugo2ens.csv')   # from BioMART
HUGO_emb_dict_sup = dict(zip(hugo_ens_sup['hgnc_symbol'],hugo_ens_sup['ensembl_gene_id']))
HUGO_emb_dict.update(HUGO_emb_dict_sup)
all_ensembl = HUGO_emb_dict.values()

gene_info_lm = gene_info['ensembl_id'][gene_info['feature_space'] == "landmark"].values
gene_info = gene_info['ensembl_id']
core_cellines = ['A375','A549','HCC515','HEPG2','HT29','MCF7','PC3','HA1E','VCAP'] # 9 core cell lines
# corresponding to tissues of: skin, lung, lung, liver, large_intestine(colon), breast, prostate  ,kidney , prostate


# 1. only parse specific interested cells

'''for interest_cell in core_cellines:
    data_df = pd.read_csv("../../CMap/Format_data/ko/ko_{}_lm.csv".format(interest_cell))
    if interest_cell == 'A375':
        data_all = data_df
    else:
        data_all = pd.concat([data_df,data_all])'''


# 2. parse data across all cells
data_all = pd.read_csv("../../CMap/Format_data/ko/ko_all_lm.csv")

data_all = pd.DataFrame(data_all)
gene_all_data = data_all['cid'].unique()
n_ge = len(gene_all_data)

print("Number of all genes in gene_info: {}".format(len(gene_info)))
print("Number of genes in data:{}".format(len(gene_all_data)))

print(data_all.shape)
count_nonexist = 0
count_noensembl = 0
count_lm = 0
flag = 0
col_name = []

maxCol=lambda x: max(x.min(), x.max(), key=abs)

for i in range(0,n_ge):
    gene_name = gene_all_data[i]
    data_i = data_all.loc[data_all['cid'] == gene_name].drop('cid',axis=1)
    # print(data_i)

    '''if data_i.empty :
        count_nonexist += 1
        continue'''
    if gene_name not in all_ensembl:
        count_noensembl += 1
        continue
    if gene_name in gene_info_lm:
        count_lm += 1


    col_name.append(gene_name)

    # print(type(data_i.iloc[2]['ENSG00000090861']))  # THIS IS FLOAT
    data_i = data_i.apply(maxCol,axis=0)

    if flag == 0:
        data_ge_maxpool = data_i
        flag = 1
    else:
        data_ge_maxpool = pd.concat([data_ge_maxpool,data_i],axis=1)

    if i % 1000 == 0:
        print(i/n_ge * 100)

data_ge_maxpool.columns = col_name
data_ge_maxpool = pd.DataFrame(data_ge_maxpool).T


# data_ge_maxpool = data_ge_maxpool.set_index('cid')
data_ge_maxpool.index.names = [None]

print(data_ge_maxpool.shape)
# print("Non_exist:{}".format(count_nonexist))
print("Non_ensembl:{}".format(count_noensembl))
print("Landmark gene:{}".format(count_lm))

'''
# diag matrix for landmark genes
lm_gene = pd.Series(list(data_all.columns)[1:])
diag = np.eye(len(lm_gene))
diag = pd.DataFrame(diag)
diag.columns = lm_gene
diag.insert(loc=0, column='gene', value= lm_gene)
diag.set_index('gene', inplace=True)
data_ge_maxpool = pd.concat([data_ge_maxpool,diag])'''

data_ge_maxpool.to_csv("../../CMap/Format_data/Maxpool/ko_all.csv",index_label=['gene'])


