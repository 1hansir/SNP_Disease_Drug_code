import pandas as pd

cp_info_d = pd.read_csv("../../CMap/Datasets/compoundinfo_beta.txt", sep="\t", low_memory=False, usecols=['pert_id', 'canonical_smiles', 'inchi_key'])
print(cp_info_d.columns)
n_ori = len(cp_info_d)
print(n_ori)

map_toDB = pd.read_csv("../../mapping_files/structure links.csv", low_memory=False, usecols=['DrugBank ID', 'InChIKey'])
map_toVA = pd.read_csv("../../mapping_files/Drugbank_rxNorm_mapping.csv", low_memory=False, usecols=['DrugBank_ID', 'RxNorm_name'])
# print(map_toDB)
# print(cp_info.loc[cp_info.duplicated()])
cp_info = cp_info_d.drop_duplicates()                                    # Duplicated compounds
print('Duplicated compounds:{}'.format(len(cp_info)-n_ori))
n_ori = len(cp_info)
cp_info = cp_info[~ cp_info['inchi_key'].isin(['restricted'])]           # restricted drugs
print('Restricted compounds:{}'.format(len(cp_info)-n_ori))
n_ori = len(cp_info)
cp_info = cp_info[~ cp_info['inchi_key'].isna()]                         # compounds that dont has corresponding INCHIKEY
print('NaN compounds:{}'.format(len(cp_info)-n_ori))
cp_info = cp_info.reset_index(drop=True)

n_cp = len(cp_info)

map_toDB_inchik = list(map_toDB['InChIKey'])
map_toDB_dbid = list(map_toDB['DrugBank ID'])

map_toVA_dbid = list(map_toVA['DrugBank_ID'])
map_toVA_RX = list(map_toVA['RxNorm_name'])

count_toDB = 0
count_toRx = 0
map_pairs = set()

for i in range(0,n_cp):
    DB_ID = ''
    Rx_name = ''
    inchi_k = cp_info.loc[i,'inchi_key']

    if inchi_k in map_toDB_inchik:
        count_toDB += 1
        k_index = map_toDB_inchik.index(inchi_k)
        DB_ID = map_toDB_dbid[k_index]

        if DB_ID in map_toVA_dbid:
            db_index = map_toVA_dbid.index(DB_ID)
            Rx_name = map_toVA_RX[db_index]
            count_toRx += 1

        map_pairs.add((cp_info.loc[i,'pert_id'],cp_info.loc[i,'canonical_smiles'],inchi_k,DB_ID,Rx_name))

map_pairs = pd.DataFrame(map_pairs)
map_pairs.to_csv("../mapping_files/pert_id2RxNorm.csv",header=['pert_id','Smiles','inchikey','DBID','RxNorm'],sep='\t')

print("# of Cleaned compounds in Cmap:{}".format(n_cp))
print("Perturbagen that can be mapped to DrugBank ID:{}".format(count_toDB))
print("Perturbagen that can be mapped to RxNorm codes:{}".format(count_toRx))





