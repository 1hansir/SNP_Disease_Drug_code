import pandas as pd

mdir = 'D:/Pycharm_projects/SNP_Disease_Drug/'
sig_info = pd.read_csv("../../CMap/Datasets/siginfo_beta.txt", sep="\t", low_memory=False)

for interest_pert in ['trt_oe','trt_xpr']:
    pert_ids = sig_info["cmap_name"][sig_info["pert_type"] == interest_pert]

    pert_ids.to_csv(mdir + 'mapping_files/HUGOSymbols_enid/HUGO_{}.csv'.format(interest_pert))