import pandas as pd
import numpy as np

mdir = 'D:/Pycharm_projects/SNP_Disease_Drug/'
gtex_filenames = ['GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Colon_Sigmoid.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Kidney_Cortex.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Liver.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Lung.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Prostate.signifpairs.txt',
                 'GTEx/GTEx_Analysis_v8_eQTL/Skin_Not_Sun_Exposed_Suprapubic.signifpairs.txt',
                  'GTEx/GTEx_Analysis_v8_eQTL/Ovary.signifpairs.txt',
                  'GTEx/GTEx_Analysis_v8_eQTL/Stomach.signifpairs.txt']

# read gtex data from different tissues:
gtex_alldata = []
for tissue_file in gtex_filenames:
    gtex_data_t = pd.read_csv(mdir + tissue_file, sep='\t')
    if tissue_file == 'GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt':
        gtex_alldata = gtex_data_t
    else:
        gtex_alldata = pd.concat([gtex_alldata, gtex_data_t])

all_SNPs = gtex_alldata['variant_id'].unique()
n_gtex = len(all_SNPs)
print("# of SNPs(unique) in GTEx file:{}".format(n_gtex))

# gtex_2_dbsnp_filename = 'mapping_files/GTExID_2_dbSNPID.lookup_table.txt'
# map_todbsnp = pd.read_csv(mdir+gtex_2_dbsnp_filename,sep='\t', low_memory=False, usecols=['variant_id', 'rs_id_dbSNP151_GRCh38p7'])
# n_map_ori = len(map_todbsnp)
# map_todbsnp = map_todbsnp[~ map_todbsnp['rs_id_dbSNP151_GRCh38p7'].isin(['.'])]   # some mapping pairs are without rsID
# print(len(map_todbsnp) - n_map_ori)
# print(len(map_todbsnp))
# map_todbsnp = map_todbsnp.reset_index(drop=True)
# map_todbsnp.to_csv(mdir+'mapping_files/GTExID_2_dbSNPID_CLEAN.lookup_table.txt',sep='\t')

gtex_2_dbsnp_filename = 'mapping_files/GTExID_2_dbSNPID_CLEAN.lookup_table.txt'
map_todbsnp = pd.read_csv(mdir+gtex_2_dbsnp_filename,sep='\t', low_memory=False, usecols=['variant_id', 'rs_id_dbSNP151_GRCh38p7'])

# Get the SNPs in GWAS catalog
gwas_filename = mdir + 'mapping_files/gwas_catalog_v1.0-associations_e109_r2023-04-24.tsv'
gwas_catalog = pd.read_csv(gwas_filename, sep = '\t', low_memory=False)
gwas_SNPs = gwas_catalog['SNPS'].unique()
print("GWAS catalog size: {}".format(len(gwas_SNPs)))


dbsnp_id_ls = list(map_todbsnp['rs_id_dbSNP151_GRCh38p7'])
variant_id_ls = list(map_todbsnp['variant_id'])

count_toGWAS = 0
count_todbSNPs = 0
sig_snp = set()


snp_db_dic = dict(zip(map_todbsnp['variant_id'],map_todbsnp['rs_id_dbSNP151_GRCh38p7']))
snp_gwas_dic = dict(zip(gwas_SNPs,gwas_SNPs))

for i in range(0,n_gtex):
    variant_id = all_SNPs[i]

    try:
        snp_id = snp_db_dic[variant_id]
        count_todbSNPs += 1

        try:
            snp_id_gwas = snp_gwas_dic[snp_id]
            count_toGWAS += 1
            sig_snp.add((variant_id, snp_id_gwas, 1))

        except KeyError:
            sig_snp.add((variant_id, snp_id, 0))

    except KeyError:
        pass

    if i % 10000 == 0:   # Will throw out division error (denominator == 0)
        print("{}%".format(100*i/n_gtex))
        print('Hit rate of VariantID to dbSNPID:{}'.format(count_todbSNPs / (i+1)))
        print('Hit rate of belong to GWAS catalog:{}'.format(count_toGWAS / (i+1)))
        print('-------------------------------------------')



print('Hit rate of VariantID to dbSNPID:{}'.format(count_todbSNPs/n_gtex))
print('Hit rate of belong to GWAS catalog:{}'.format(count_toGWAS/n_gtex))

sig_snp = pd.DataFrame(sig_snp)
sig_snp.columns = ['variant_id','SNP_id','GWAS_label']
sig_snp.to_csv(mdir+'mapping_files/SNP_INDEX/snpid_file.csv',sep='\t',header=['variant_id','SNP_id','GWAS_label'],index_label='SNP_index')

sig_snp_ingwas = sig_snp[sig_snp['GWAS_label'] == 1]
sig_snp_ingwas.to_csv(mdir+'mapping_files/SNP_INDEX/snpid_GWAS_file.csv',sep='\t',columns=['variant_id','SNP_id'],
                      header=['variant_id','SNP_id'],index_label='SNP_index')

