# remove all variables
print(getwd())
rm(list = ls())

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

Label_ = data.frame(dataSet$sample.group,
                    get_true_num_label(dataSet$sample.group))

true = get_true_num_label(dataSet$sample.group)
true_char = dataSet$sample.group

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
# test different methods
res_out = list(); t=1
ks = c(500, 1000, 2000, 3000, 4000, 5000)
avgNMI_Mat = matrix(0, 4, length(ks))
colnames(avgNMI_Mat) = paste("k =",ks)
row.names(avgNMI_Mat) = c("nmf_l20", "ONMF_l20", "ONMF_l20_rho","Kmeans_l0")

for(k in ks){
  ResMat = matrix(0, 4, 10)  
  for(j in 1:10){
    print(j)
    # initialize 
    W0H0 = get_unif_init(X, ncl, nstart0=j)
    Winit = W0H0$W; Hinit = W0H0$H
    
    source('../My_code_myFuns/NMF_BCD.R')
    out0 = NMF(X, Winit, Hinit, tol=1e-3, maxit=500)
    
    # NMF_l20 algorithm
    source('../My_code_myFuns/SNMF_L20_new.R')
    out7 = SNMF_L20(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
    NMI7 = get_NMI(true_char, apply(out7$H,2,which.max))
    
    source('../My_code_myFuns/SONMF_L20_new.R')
    out2 = SONMF_L20(X, Winit, Hinit, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
    NMI2 = get_NMI(true_char, get_optW_optH(out2)$preLab)
    
    source('../My_code_myFuns/SONMF_L20_FixRho_new.R')
    out3 = SONMF_L20_FixRho(X, Winit, Hinit, k, rho=1, palm_tol=1e-3, palm_maxiter)
    NMI3 = get_NMI(true_char, apply(out3$H,2,which.max))
    
    # sparse K-means
    source('../My_code_myFuns/sparse_kmeans.R')
    set.seed(j)
    out10 = L0_kmeans(t(X), ncl, k, nstart = 1)
    # NMI10 = external_validation(true, out10$cluster, method= "nmi")
    NMI10 = get_NMI(true_char, out10$cluster)
    
    ResMat[,j] = c(NMI7, NMI2, NMI3, NMI10)
  }
  row.names(ResMat) = c("nmf_l20", "ONMF_l20", "ONMF_l20_rho","Kmeans_l0")
  res_out[[t]] = ResMat
  avgNMI_Mat[,t] = apply(ResMat,1,mean)
  t = t + 1
  
  print(avgNMI_Mat)
}

# NMIs = c()
# for(i in 1:length(res_out)){
#   NMIs[i] = mean(res_out[[i]]["nmi",])
# }

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# get LaTex Table 
Tab1 = format(round(avgNMI_Mat*100,2), nsmall = 2)
library(xtable)
xtable(Tab1)

time = proc.time() - ptm; print(time)
save(res_out, avgNMI_Mat, ks, Tab1, time, file = "F12_pollen.RData")




