install.packages("anndata")
install.packages(c("SeuratDisk", "sceasy"))
library(anndata)


adata <- read_h5ad("/Users/mahimasrivastava/Hd5 files/GSM4186971.h5ad")
print(adata)
