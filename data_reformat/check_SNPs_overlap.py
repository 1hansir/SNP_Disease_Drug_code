import pandas as pd
import numpy as np

mdir = 'D:/Pycharm_projects/SNP_Disease_Drug/'
gtex_filename = mdir + 'GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt'

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

gtex_alldata = pd.read_csv(gtex_filename,sep='\t')
all_SNPs = gtex_alldata['variant_id'].unique()
n_gtex = len(all_SNPs)
print("# of SNPs(unique) in GTEx file:{}".format(n_gtex))

dbsnp_id_ls = list(map_todbsnp['rs_id_dbSNP151_GRCh38p7'])
variant_id_ls = list(map_todbsnp['variant_id'])

count_toGWAS = 0
count_todbSNPs = 0
sig_snp = set()
gwas_mapped = []
gwas_unmapped = []
gtex_unmapped = []
for i in range(0,n_gtex):
    variant_id = all_SNPs[i]

    try:
        index_dbsnp = variant_id_ls.index(variant_id)
        count_todbSNPs += 1
        snp_id = dbsnp_id_ls[index_dbsnp]

        if snp_id in gwas_SNPs:
            count_toGWAS += 1
            sig_snp.add((variant_id, snp_id))
            gwas_mapped.append(snp_id)
        else:
            gtex_unmapped.append(snp_id)

    except ValueError:
        pass

    if i % 10000 == 0:   # Will throw out division error (denominator == 0)
        print("{}%".format(100*i/n_gtex))
        print('Hit rate of VariantID to dbSNPID:{}'.format(count_todbSNPs / (i+1)))
        print('Hit rate of belong to GWAS catalog:{}'.format(count_toGWAS / (i+1)))
        print('-------------------------------------------')



print('Hit rate of VariantID to dbSNPID:{}'.format(count_todbSNPs/n_gtex))
print('Hit rate of belong to GWAS catalog:{}'.format(count_toGWAS/n_gtex))

sig_snp = pd.DataFrame(sig_snp)
sig_snp.to_csv(mdir+'mapping_files/sig_snpid_inGWAS.csv',sep='\t',header=['variant_id','SNP_id'])

gwas_mapped = pd.DataFrame(gwas_mapped)
gwas_mapped.to_csv(mdir+'mapping_files/gwas_mapped.csv',sep='\t')
gtex_unmapped = pd.DataFrame(gtex_unmapped)
gtex_unmapped.to_csv(mdir+'mapping_files/gtex_unmapped.csv',sep='\t')

# gwas_unmapped = list(set(gwas_SNPs)-set(gwas_mapped))
