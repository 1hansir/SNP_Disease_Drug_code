from cmapPy.pandasGEXpress.parse import parse
import pandas as pd

# cell_info = pd.read_csv("../../CMap/Datasets/cellinfo_beta.txt", sep="\t", low_memory=False)
# sig_info = pd.read_csv("../../CMap/Datasets/siginfo_beta.txt", sep="\t", low_memory=False)
# gene_info = pd.read_csv("../../CMap/Datasets/geneinfo_beta.txt", sep="\t", dtype=str)
cp_info_d = pd.read_csv("../../CMap/Datasets/compoundinfo_beta.txt", sep="\t", low_memory=False, usecols=['pert_id'])
cp_info_d = cp_info_d['pert_id']
n_cp = len(cp_info_d)
core_cellines = ['A375','A549','HCC515','HEPG2','HT29','MCF7','PC3','HA1E','VCAP'] # 9 core cell lines
# corresponding to tissues of: skin, lung, lung, liver, large_intestine(colon), breast, prostate  ,kidney , prostate
col_name = []


for interest_cell in core_cellines:   # concatenate all the data
    data_df = pd.read_csv("../../CMap/Format_data/cp/cp_{}_lm.csv".format(interest_cell))
    if interest_cell == 'A375':
        data_all = data_df
    else:
        data_all = pd.concat([data_df,data_all])

data_all = pd.DataFrame(data_all)

print(data_all.shape)
count_nonexist = 0
flag = 0
data_cp_maxpool = []
maxCol=lambda x: max(x.min(), x.max(), key=abs)


# max pooling of every compound's data across different cell lines
for i in range(0,n_cp):
    cp_name = cp_info_d[i]
    data_i = data_all.loc[data_all['cid'] == cp_name].drop('cid',axis=1)
    # print(data_i)

    if data_i.empty:
        count_nonexist += 1
        continue

    col_name.append(cp_name)
    data_i = data_i.apply(maxCol,axis=0)

    if flag == 0:
        data_cp_maxpool = data_i
        flag = 1
    else:
        data_cp_maxpool = pd.concat([data_cp_maxpool,data_i],axis=1)

    if i % 1000 == 0:
        print(i/n_cp * 100)


data_cp_maxpool.columns = col_name
data_cp_maxpool = pd.DataFrame(data_cp_maxpool).T


print(data_cp_maxpool.shape)
print(count_nonexist)
data_cp_maxpool.to_csv("../../CMap/Format_data/Maxpool/cp.csv", index_label=['compound'])