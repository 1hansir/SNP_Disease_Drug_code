import pandas as pd

cp_info_d = pd.read_csv("../../CMap/Datasets/compoundinfo_beta.txt", sep="\t", low_memory=False, usecols=['pert_id', 'canonical_smiles', 'inchi_key'])
print(cp_info_d.columns)
n_ori = len(cp_info_d)
print(n_ori)

map_toDB = pd.read_csv("../../mapping_files/Drug_Bank/drugbank vocabulary.csv", low_memory=False,sep=',', usecols=['DrugBank ID', 'Standard InChI Key'])
# map_toDB_2 = pd.read_csv("../../mapping_files/Drug_Bank/structure links.csv", low_memory=False, usecols=['DrugBank ID', 'InChIKey'])
map_toVA = pd.read_csv("../../mapping_files/DB_RN_Map/RxNorm_DB.csv", low_memory=False, usecols=['RXCUI', 'CODE'])
map_toDB.columns = ['DrugBank ID', 'InChIKey']
map_toVA.columns = ['RxNorm_name','DrugBank_ID']


VA_filename = 'embeddings/VA_embed_500_20220607_items.csv'
VA_item = pd.read_csv('../../embeddings/VA_embed_500_20220607_items.csv',usecols=['codes'])
VA_RX = VA_item[VA_item['codes'].str.startswith('RXNORM:',na=False)]
VA_RX = VA_RX.reset_index(drop=True)

# sub string from RXNORM: to pure code
for i in range(0,len(VA_RX)):
    rx_code = VA_RX.loc[i].str
    code = rx_code[7:]
    VA_RX.loc[i] = code

print('RXNORM codes in VA EHR space:{}'.format(len(VA_RX)))

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

map_toDB = map_toDB[~map_toDB['InChIKey'].isna()]
map_toDB_inchik = list(map_toDB['InChIKey'])
map_toDB_dbid = list(map_toDB['DrugBank ID'])

map_toVA_dbid = list(map_toVA['DrugBank_ID'])
map_toVA_RX = list(map_toVA['RxNorm_name'])

map_toVA = list(VA_RX['codes'])

count_toDB = 0
count_toRx = 0
count_toVA = 0
map_pairs = set()
map_pairs_toVA = set()

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

            if str(Rx_name) in map_toVA:
                # RX_index = map_toVA.index(Rx_name)
                # Rx_name = map_toVA[RX_index]
                count_toVA += 1
                map_pairs_toVA.add((cp_info.loc[i, 'pert_id'], Rx_name, DB_ID,  cp_info.loc[i, 'canonical_smiles'], inchi_k))

        map_pairs.add((cp_info.loc[i,'pert_id'],Rx_name,DB_ID,cp_info.loc[i,'canonical_smiles'],inchi_k))

map_pairs = pd.DataFrame(map_pairs)
map_pairs.to_csv("../../mapping_files/pert_id2RxNorm.csv",header=['pert_id','RxNorm','DBID','Smiles','inchikey'],sep='\t')
map_pairs_toVA = pd.DataFrame(map_pairs_toVA)
map_pairs_toVA.to_csv("../../mapping_files/pert_id2VA.csv",header=['pert_id','RxNorm','DBID','Smiles','inchikey'],sep='\t')

print("# of Cleaned compounds in Cmap:{}".format(n_cp))
print("Perturbagen that can be mapped to DrugBank ID:{}".format(count_toDB))
print("Perturbagen that can be mapped to RxNorm codes through DB:{}".format(count_toRx))
print("Perturbagen that can be mapped to VA Space:{}".format(count_toVA))





