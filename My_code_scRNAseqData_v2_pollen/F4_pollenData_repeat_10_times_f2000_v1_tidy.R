# ------------------------------------------------------------------------------
# remove all variables
rm(list = ls())
# ------------------------------------------------------------------------------
# R packages
library("ggplot2","ClusterR","NMI")

# import helper functions
source('../My_code_myFuns/Add_Function_FULL.R')
# ------------------------------------------------------------------------------
# get time
ptm <- proc.time()
# ------------------------------------------------------------------------------
# get data
# 1_pollenData
# Samples: 301
# Genes: 8747
# Classes: 4
load("../My_data/1_pollenData/pollenData.RData")
dataSet = pollenData

Exp = dataSet$Exp
X = log2(Exp+1)
true = dataSet$sample.group
dat_label = data.frame(dataSet$sample.group,
                    get_true_num_label(dataSet$sample.group))
ncl = length(unique(true))

A = as.matrix(dataSet$ppi_Net)
# ------------------------------------------------------------------------------
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
ResMat = matrix(0, 11, repeatNum)  
for(j in 1:repeatNum){
  # initialize 
  print(j)
  W0H0 = get_unif_init(X, ncl, nstart0=j)
  Winit = W0H0$W; Hinit = W0H0$H
  
  source('../My_code_myFuns/NMF_BCD.R')
  out0 = NMF(X, Winit, Hinit, tol=1e-3, maxit=500)
  
  source('../My_code_myFuns/ONMF_new.R')
  out1 = ONMF(X, Winit, Hinit, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
  
  source('../My_code_myFuns/SONMF_L20_new.R')
  out2 = SONMF_L20(X, Winit, Hinit, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
  
  source('../My_code_myFuns/SONMF_L20_FixRho_new.R')
  out3 = SONMF_L20_FixRho(X, Winit, Hinit, k, rho=1, palm_tol=1e-3, palm_maxiter)
  
  source('../My_code_myFuns/SONMF_CL0_new.R')
  out4 = SONMF_CL0(X, out0$W, out0$H, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
  
  source('../My_code_myFuns/SONMF_L0_new.R')
  out5 = SONMF_L0(X, out0$W, out0$H, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
  
  source('../My_code_myFuns/SNMF_L0_new.R')
  out6 = SNMF_L0(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
  
  source('../My_code_myFuns/SNMF_L20_new.R')
  out7 = SNMF_L20(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
  
  source('../My_code_myFuns/SNMF_CL0_new.R')
  out8 = SNMF_CL0(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
  
  # K-means
  set.seed(j)
  out9 = kmeans(t(X), ncl, nstart = 1)
  # out9 = kmeans(t(X), t(Winit))
  NMI9 = external_validation(true, out9$cluster, method= "nmi")
  
  # sparse K-means
  #source('../My_code_myFuns/L0_kmeans.R')
  source('../My_code_myFuns/sparse_kmeans.R')
  set.seed(j)
  out10 = L0_kmeans(t(X), ncl, k, nstart = 1)
  
  # save the result
  outs = list(out0,out1,out2,out3,out4,out5,out6,out7,out8,out9,out10)
  names(outs) = c("NMF", "ONMF","ONMF_L20", "ONMF_L20_Rho","ONMF_CL0",
                  "ONMF_L0","NMF_L0","NMF_L20","NMF_CL0","Kmeans","L0_Kmeans")
  outs_list[[j]] = list(out0,out1,out2,out3,out4,out5,out6,out7,out8,out9,out10)
}

# ------------------------------------------------------------------------------
# get time
time = proc.time() - ptm; print(time)

fileName = paste("F4_res_pollenData_repeat_10_times_f2000", ".RData", sep = "")  
save(outs_list, dat_label, file = fileName)
# ------------------------------------------------------------------------------