import pandas as pd
import numpy as np

cell_info = pd.read_csv("../../CMap/Datasets/cellinfo_beta.txt", sep="\t", low_memory=False)

print(cell_info.columns)
print(cell_info['cell_iname'][cell_info['cell_lineage']=='breast'])
print(cell_info['cell_iname'][cell_info['cell_lineage']=='large_intestine'])
print(cell_info['cell_iname'][cell_info['cell_lineage']=='kidney'])
print(cell_info['cell_iname'][cell_info['cell_lineage']=='liver'])
print(cell_info['cell_iname'][cell_info['cell_lineage']=='lung'])
print(cell_info['cell_iname'][cell_info['cell_lineage']=='prostate'])
print(cell_info['cell_iname'][cell_info['cell_lineage']=='skin'])
print(cell_info['cell_iname'][cell_info['cell_lineage']=='ovary'])
print(cell_info['cell_iname'][cell_info['cell_lineage']=='stomach'])


data_all = pd.read_csv("../../CMap/Format_data/oe/oe_all_lm.csv")
lm_gene = pd.Series(list(data_all.columns)[1:])
diag = np.eye(len(lm_gene))
diag = pd.DataFrame(diag)
diag.columns = lm_gene
diag.insert(loc=0, column='index', value= lm_gene)
diag.set_index('index', inplace=True)
print(diag)

