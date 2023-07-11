library(data.table)
library(dplyr)

mdir <- 'D:/Pycharm_projects/SNP_Disease_Drug/'

cmap_filename <- paste0(mdir, '../../CMap/Format_data/cp_MCF7_lm.csv')
cmap_alldata <- fread(cmap_filename)
lm_gene_ls <- unlist(colnames(cmap_alldata))
lm_gene_ls <- lm_gene_ls[2:length(lm_gene_ls)]

gtex_filename <- paste0(mdir, '../../GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.signifpairs.txt')
# gtex_filename <- paste0(mdir, 'GTEx/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.v8.egenes.txt')
gtex_alldata <- fread(gtex_filename)

# Method1: first select p-value, then select Gene ID
gtex_alldata <- gtex_alldata[abs(gtex_alldata$pval_nominal)<1e-30]

n_gtex <- length(gtex_alldata$variant_id)

all_slope <- gtex_alldata$slope
all_SNPs <- gtex_alldata$variant_id
all_GENEs <- gtex_alldata$gene_id

print(length(all_SNPs))
uni_SNPs <- unique(all_SNPs)
print(length(uni_SNPs))

all_GENEs <- unlist(lapply(all_GENEs,FUN=function(x){unlist(strsplit(x,"\\."))[1]}))
print(length(all_GENEs))
uni_GENEs <- unique(all_GENEs)
print(length(uni_GENEs))
# EXAMINE the number of Genes and SNPs


# represent genes with correlation with different SNPs:
all_df <- data.frame(matrix(data = 0, nrow = length(lm_gene_ls),ncol = length(uni_SNPs)),row.names = lm_gene_ls)
names(all_df) <- c(uni_SNPs)

count_hit <- 0
for (i in 1:n_gtex){
  if (all_GENEs[i] %in% lm_gene_ls){
      all_df[all_GENEs[i],all_SNPs[i]] <- all_slope[i]
      count_hit <- count_hit + 1
  }
}

# Method 2: First select gene, then select p-value
all_GENEs <- gtex_alldata$gene_id
all_GENEs <- unlist(lapply(all_GENEs,FUN=function(x){unlist(strsplit(x,"\\."))[1]}))
gtex_alldata$gene_id <- data.frame(all_GENEs)

gtex_alldata <- gtex_alldata[gtex_alldata$gene_id %in% lm_gene_ls]
gtex_alldata <- gtex_alldata[abs(gtex_alldata$pval_nominal)<1e-10]

n_gtex <- length(gtex_alldata$variant_id)

all_slope <- gtex_alldata$slope
all_SNPs <- gtex_alldata$variant_id
all_GENEs <- gtex_alldata$gene_id

print(length(all_SNPs))
uni_SNPs <- unique(all_SNPs)
print(length(uni_SNPs))

print(length(all_GENEs))
uni_GENEs <- unique(all_GENEs)
print(length(uni_GENEs))

all_df <- data.frame(matrix(data = 0, nrow = length(lm_gene_ls),ncol = length(uni_SNPs)),row.names = lm_gene_ls)
names(all_df) <- c(uni_SNPs)

for (i in 1:n_gtex){
  all_df[all_GENEs[i],all_SNPs[i]] <- all_slope[i]
}

fwrite(all_df, paste0(mdir, '../../GTEx/Format_data/lmgene.csv'))








