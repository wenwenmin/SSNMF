# remove all variables
rm(list = ls())
# ------------------------------------------------------------------------------
# R packages
library("ggplot2")
library("ClusterR")
library("NMI")

# import helper functions
source('../My_code_myFuns/Add_Function_FULL.R')
# ------------------------------------------------------------------------------

# get time
ptm <- proc.time()

# ------------------------------------------------------------------------------
# load data
# 1_pollenData
# Samples: 301
# Genes: 8747
# Classes: 4
load("../My_data/1_pollenData/pollenData.RData")
dataSet = pollenData

Exp = dataSet$Exp
X = log2(Exp+1)
# true = dataSet$sample.group
true = get_true_num_label(dataSet$sample.group) 

Label_ = data.frame(dataSet$sample.group,
                    get_true_num_label(dataSet$sample.group))

ncl = length(unique(true))

A = as.matrix(dataSet$ppi_Net)

# ------------------------------------------------------------------------------
# parameter setting
k = 2000 
rho = .1 
gama = 1.5
palm_tol = 1e-3 
rho_tol = 1e-3
palm_maxiter = 100
rho_maxiter = 10
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# test different methods
outs_list = list()
repeatNum = 10
ResMat = matrix(0, 3, repeatNum)  
gene_num_vec = c()

for(j in 1:repeatNum){
  print(j)
  
  # # K-means
  # set.seed(j)
  # out9 = kmeans(t(X), ncl, nstart = 1)
  # NMI9 = external_validation(true, out9$cluster, method= "nmi")
  # 
  # # sparse K-means
  # source('../My_code_myFuns/sparse_kmeans.R')
  # set.seed(j)
  # out10 = L0_kmeans(t(X), ncl, k, nstart = 1)
  # NMI10 = external_validation(true, out10$cluster, method= "nmi")
  
  # sparse K-means with L1 penalty
  library("sparcl")
  set.seed(j)
  r11 = KMeansSparseCluster(t(X), ncl, wbounds=27, nstart = 5)

  gene_num_vec[j] = length(which(r11[[1]]$ws!=0))
  
  NMI =     external_validation(true, r11[[1]]$Cs, method= "nmi")
  purity =  external_validation(true, r11[[1]]$Cs, method= "purity")
  entropy = external_validation(true, r11[[1]]$Cs, method= "entropy")
  
  # save results
  ResMat[,j] = c(NMI, purity, entropy)
}

# ------------------------------------------------------------------------------
res = c(apply(ResMat,1,mean), mean(gene_num_vec))
names(res) = c("nmi", "purity", "entropy","gene_num")
print(res)
# nmi       purity      entropy     gene_num 
# 0.7809497    0.8853821    0.2145277 2039.0000000 

# ------------------------------------------------------------------------------
avg = apply(ResMat,1,mean)
sd  = apply(ResMat,1,sd)
t1 = format(round(avg*100,2), nsmall = 2)
t2 = format(round(sd*100,2), nsmall = 2)

rowRes = c(mean(gene_num_vec), paste(t1,"+",t2))
names(rowRes) = c("gene_num", "nmi", "purity", "entropy")

print(rowRes)

# ------------------------------------------------------------------------------
rm(A,dataSet,Exp)
time = proc.time() - ptm; print(time)

# ------------------------------------------------------------------------------
fileName = paste("F10_pollen_", "v1.RData", sep = "")  
save.image(fileName)




