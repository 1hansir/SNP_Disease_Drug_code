import matplotlib.pyplot as plt
import pickle

mdir = 'D:/Pycharm_projects/SNP_Disease_Drug/'

# Plot the histogram distribution of SNP-Gene pairs
with open(mdir + 'GTEx/Format_data/saved_snp_dict.pkl', 'rb') as f:
    loaded_dict = pickle.load(f)

with open(mdir + 'GTEx/Format_data/saved_gene_dict.pkl', 'rb') as f:
    pair_dic = pickle.load(f)
with open(mdir + 'GTEx/Format_data/saved_snp_dict.pkl', 'rb') as f:
    snp_pair_dic = pickle.load(f)

with open(mdir + 'GTEx/Format_data/saved_gene_num.pkl', 'rb') as f:
    pair_dic_num = pickle.load(f)
with open(mdir + 'GTEx/Format_data/saved_snp_num.pkl', 'rb') as f:
    snp_pair_dic_num = pickle.load(f)

print(len(pair_dic))
print(max(pair_dic_num.values()))


print(len(snp_pair_dic))
print(max(snp_pair_dic_num.values()))

n, bins, patches = plt.hist(pair_dic_num.values(),range=(0,1000),bins=100,log=False)
for i in range(len(n)):
    plt.text(bins[i]+(bins[1]-bins[0])/2, n[i]*1.01, int(n[i]), ha='center', va= 'bottom')
plt.show()

n, bins, patches = plt.hist(snp_pair_dic_num.values(),bins=19,log=False)
for i in range(len(n)):
    plt.text(bins[i]+(bins[1]-bins[0])/2, n[i]*1.01, int(n[i]), ha='center', va= 'bottom')
plt.show()

print('test')

