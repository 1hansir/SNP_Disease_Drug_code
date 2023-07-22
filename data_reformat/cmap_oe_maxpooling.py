from cmapPy.pandasGEXpress.parse import parse
import pandas as pd

# cell_info = pd.read_csv("../../CMap/Datasets/cellinfo_beta.txt", sep="\t", low_memory=False)
# sig_info = pd.read_csv("../../CMap/Datasets/siginfo_beta.txt", sep="\t", low_memory=False)
gene_info = pd.read_csv("../../CMap/Datasets/geneinfo_beta.txt", sep="\t", dtype=str)
# cp_info_d = pd.read_csv("../../CMap/Datasets/compoundinfo_beta.txt", sep="\t", low_memory=False, usecols=['pert_id'])

HUGO_emb_dict = dict(zip(gene_info['gene_symbol'],gene_info['ensembl_id']))
gene_info_lm = gene_info['ensembl_id'][gene_info['feature_space'] == "landmark"]
gene_info = gene_info['ensembl_id']
n_ge = len(gene_info)
core_cellines = ['A375','A549','HCC515','HEPG2','HT29','MCF7','PC3','HA1E','VCAP'] # 9 core cell lines
# corresponding to tissues of: skin, lung, lung, liver, large_intestine(colon), breast, prostate  ,kidney , prostate



for interest_cell in core_cellines:
    data_df = pd.read_csv("../../CMap/Format_data/oe/oe_{}_lm.csv".format(interest_cell))
    if interest_cell == 'A375':
        data_all = data_df
    else:
        data_all = pd.concat([data_df,data_all])

data_all = pd.DataFrame(data_all)
gene_all_data = data_all['cid'].unique()

print("Number of all genes in gene_info: {}".format(n_ge))
print("Number of genes in data:{}".format(len(gene_all_data)))

print(data_all.shape)
count_nonexist = 0
count_noensembl = 0
count_lm = 0
flag = 0
col_name = []

maxCol=lambda x: max(x.min(), x.max(), key=abs)

for i in range(0,n_ge):
    gene_name = gene_info[i]
    data_i = data_all.loc[data_all['cid'] == gene_name].drop('cid',axis=1)
    # print(data_i)

    if data_i.empty :
        count_nonexist += 1
        continue
    if gene_name not in HUGO_emb_dict.values():
        count_noensembl += 1
        continue
    if gene_name in gene_info_lm:
        count_lm += 1


    col_name.append(gene_name)

    # print(type(data_i.iloc[2]['ENSG00000090861']))  # THIS IS NUMPY.FLOAT64
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

print(data_ge_maxpool.shape)
print("Non_exist:{}".format(count_nonexist))
print("Landmark gene:{}".format(count_lm))  # = 0, make sense because pertubation should be different from what we use as expression vectors
data_ge_maxpool.to_csv("../../CMap/Format_data/Maxpool/oe.csv",index_label=['gene'])


