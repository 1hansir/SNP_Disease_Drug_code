import pkg_resources
import sys
# Print version of cmapPy being used in current conda environment
print(pkg_resources.get_distribution("cmapPy").version)

import pandas as pd

cell_info = pd.read_csv("../../CMap/Datasets/cellinfo_beta.txt", sep="\t", low_memory=False)
core_cellines = ['A375','A549','HCC515','HEPG2','HT29','MCF7','PC3','HA1E','VCAP']  # 9 core cell lines
# corresponding to tissues of: skin, lung, lung, liver, large_intestine(colon), breast, prostate  ,kidney , prostate
# (there's also endometrium, nervous_system, overy, stomach that are also in GTEx's studied tissues)

sig_info = pd.read_csv("../../CMap/Datasets/siginfo_beta.txt", sep="\t", low_memory=False)
all_cells = sig_info['cell_iname'].unique()
print(all_cells)
print(sig_info.columns)
cp_ids = sig_info["sig_id"][sig_info["pert_type"] == "trt_cp"]
cell_ids = sig_info["sig_id"][sig_info["cell_iname"] == "MCF7"]

cp_MCF7_ids = list(set(cp_ids) & set(cell_ids))
print(len(cp_MCF7_ids))


gene_info = pd.read_csv("../../CMap/Datasets/geneinfo_beta.txt", sep="\t", dtype=str)
landmark_ids = gene_info['gene_id'][gene_info['feature_space'] == "landmark"]
gene_emsembl_id = gene_info['ensembl_id'][gene_info['feature_space'] == "landmark"]
print(gene_info.columns)
print(len(landmark_ids))


from cmapPy.pandasGEXpress.parse import parse

my_col_metadata = parse("../../CMap/Datasets/level5_beta_trt_cp_n720216x12328.gctx", col_meta_only=True)
print(my_col_metadata.shape)

my_row_metadata = parse("../../CMap/Datasets/level5_beta_trt_cp_n720216x12328.gctx", row_meta_only=True)
print(my_row_metadata.shape)

only_gtxtoo = parse("../../CMap/Datasets/level5_beta_trt_cp_n720216x12328.gctx", cid=cp_MCF7_ids, rid=landmark_ids)
print(only_gtxtoo.data_df.shape)
only_gtxtoo = only_gtxtoo.data_df.T
only_gtxtoo.columns = gene_emsembl_id

only_gtxtoo.to_csv("../CMap/Format_data/cp_MCF7_lm.csv")

