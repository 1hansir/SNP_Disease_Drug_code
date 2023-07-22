import pandas as pd

mdir = 'D:/Pycharm_projects/SNP_Disease_Drug/'
rxnorm_filename =  'mapping_files/RxNorm/rrf/RXNCONSO.RRF'

rxnorm = pd.read_csv(mdir+rxnorm_filename,sep='|', low_memory=False, header=None)

rxnorm.columns = ['RXCUI','LAT','TS','LUI','STT','SUI','ISPREF','RXAUI','SAUI','SCUI','SDUI','SAB','TTY','CODE','STR','SRL','SUPPRESS','CVF','p',]

rxnorm_source = rxnorm[['RXCUI','SAB','CODE']]
rxnorm_db = rxnorm_source.loc[rxnorm_source['SAB'] == 'DRUGBANK']
rxnorm_db = rxnorm_db.reset_index(drop=True)[['RXCUI','CODE']]
rxnorm_db = rxnorm_db.drop_duplicates()

rxnorm_db.to_csv(mdir+'mapping_files/RxNorm/RxNorm_DB.csv')
print(rxnorm_db)