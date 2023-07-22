import pkg_resources
from cmapPy.pandasGEXpress.parse import parse
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
gene_info = pd.read_csv("../../CMap/Datasets/geneinfo_beta.txt", sep="\t", dtype=str)
landmark_ids = gene_info['gene_id'][gene_info['feature_space'] == "landmark"]
gene_emsembl_id = gene_info['ensembl_id'][gene_info['feature_space'] == "landmark"]


pert_type_dict = {'trt_cp':'cp','trt_oe':'oe','trt_xpr':'ko'}
pert_cp_id_dict = dict(zip(sig_info["sig_id"],sig_info["pert_id"]))
pert_gene_id_dict = dict(zip(sig_info["sig_id"],sig_info["cmap_name"]))
HUGO_emb_dict = dict(zip(gene_info['gene_symbol'],gene_info['ensembl_id']))

for interest_cellname in core_cellines:
    for interest_pert in ['trt_cp','trt_oe','trt_xpr']:
        pert_ids = sig_info["sig_id"][sig_info["pert_type"] == interest_pert]
        cell_ids = sig_info["sig_id"][sig_info["cell_iname"] == interest_cellname]
        pert_cell_speci_ids = list(set(pert_ids) & set(cell_ids))

        if interest_pert == 'trt_cp':
            pert_names = {id: pert_cp_id_dict[id] for id in pert_cell_speci_ids}
        else:
            # pert_names = {id: pert_gene_id_dict[id] for id in pert_cell_speci_ids } # cmap name
            # HUGO_emb_ids = sig_info["sig_id"][sig_info["cmap_name"] in HUGO_emb_dict.keys()]
            # pert_cell_speci_ids = list(set(pert_ids) & set(cell_ids) & set(HUGO_emb_ids))

            pert_names = {id: HUGO_emb_dict[pert_gene_id_dict[id]] for id in pert_cell_speci_ids if
                          pert_gene_id_dict[id] in HUGO_emb_dict.keys()}
        # Dictionary of index_name mapping

        # my_col_metadata = parse("../../CMap/Datasets/level5_beta_{}.gctx".format(interest_pert), col_meta_only=True)
        # print(my_col_metadata.shape)

        # my_row_metadata = parse("../../CMap/Datasets/level5_beta_{}.gctx".format(interest_pert), row_meta_only=True)
        # print(my_row_metadata.shape)

        only_gtxtoo = parse("../../CMap/Datasets/level5_beta_{}.gctx".format(interest_pert), cid=pert_cell_speci_ids, rid=landmark_ids)
        print("The shape of pert {} at cell {} is {}".format(interest_pert,interest_cellname,only_gtxtoo.data_df.shape))
        only_gtxtoo = only_gtxtoo.data_df.T



        only_gtxtoo.columns = gene_emsembl_id     # rename the columns with gene emsemble ID
        only_gtxtoo.rename(index = pert_names, inplace=True)           # rename the row indices to perturbation ids

        only_gtxtoo.to_csv("../../CMap/Format_data/{}/{}_{}_lm.csv".format(pert_type_dict[interest_pert],pert_type_dict[interest_pert], interest_cellname))

