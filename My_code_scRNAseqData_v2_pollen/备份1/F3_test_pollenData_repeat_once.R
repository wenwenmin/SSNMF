# ------------------------------------------------------------------------------
# remove all variables
rm(list = ls())
# ------------------------------------------------------------------------------
# R packages
library("ggplot2")
library("NMI")
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
ncl = length(unique(true))

A = as.matrix(dataSet$ppi_Net)
# ------------------------------------------------------------------------------
# initialize 
W0H0 = get_unif_init(X,ncl,nstart0=3)
#W0H0 = get_kmeans_init(X,ncl,nstart0=1)
Winit = W0H0$W; Hinit = W0H0$H

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
source('../My_code_myFuns/NMF_BCD.R')
out0 = NMF(X, Winit, Hinit, tol=1e-3, maxit=500)
NMI0 = get_NMI_Results(list(list(out0)))

source('../My_code_myFuns/ONMF_new.R')
out1 = ONMF(X, Winit, Hinit, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
NMI1 = get_NMI_Results(list(out1))

source('../My_code_myFuns/SONMF_L20_new.R')
out2 = SONMF_L20(X, Winit, Hinit, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
NMI2 = get_NMI_Results(list(out2))

source('../My_code_myFuns/SONMF_L20_FixRho_new.R')
ou3 = SONMF_L20_FixRho(X, Winit, Hinit, k, rho=1, palm_tol=1e-3, palm_maxiter)
NMI3 = get_NMI_Results(list(list(ou3)))

source('../My_code_myFuns/SONMF_CL0_new.R')
ou4 = SONMF_CL0(X, out0$W, out0$H, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
NMI4 = get_NMI_Results(list(ou4))

source('../My_code_myFuns/SONMF_L0_new.R')
ou5 = SONMF_L0(X, out0$W, out0$H, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
NMI5 = get_NMI_Results(list(ou5))

source('../My_code_myFuns/SNMF_L0_new.R')
out6 = SNMF_L0(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
NMI6 = get_NMI_Results(list(list(out6)))

source('../My_code_myFuns/SNMF_L20_new.R')
out7 = SNMF_L20(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
NMI7 = get_NMI_Results(list(list(out7)))

source('../My_code_myFuns/SNMF_CL0_new.R')
out8 = SNMF_CL0(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
NMI8 = get_NMI_Results(list(list(out8)))

# K-means
out9 = kmeans(t(X), ncl, nstart = 1)
NMI9 = get_NMI(true,out9$cluster)

source('../My_code_myFuns/L0_kmeans.R')
out10 = L0_kmeans(t(X), ncl, k, nstart0=1)
NMI10 = get_NMI(true,out10$cluster)

# ------------------------------------------------------------------------------
# 结果展示
NMIs_diffMet = c(NMI0, NMI1, NMI2, NMI3, NMI4, NMI5, NMI6, NMI7, NMI8, NMI9, NMI10)
names(NMIs_diffMet) = c("NMF", 
                        "ONMF",
                        "ONMF_L20", 
                        "ONMF_L20_Rho",
                        "ONMF_CL0",
                        "ONMF_L0",
                        "NMF_L0",
                        "NMF_L20",
                        "NMF_CL0",
                        "Kmeans","L0_Kmeans")

print(as.matrix(NMIs_diffMet))
# NMF          0.8313297
# ONMF         0.8515611
# ONMF_L20     0.8333502
# ONMF_L20_Rho 0.8678727
# ONMF_CL0     0.8151869
# ONMF_L0      0.5284175
# NMF_L0       0.5514304
# NMF_L20      0.8678727
# NMF_CL0      0.8167976
# Kmeans       0.8251795
# L0_Kmeans    0.6477611

# ------------------------------------------------------------------------------
# 生物信息分析[待完成]
# get Network enrichment result
# ResGeneEdges = get_geneEdge_Results(Results)

# (4) get time
time = proc.time() - ptm; print(time)
rm(dataSet)
fileName = paste("F3_res_pollenData", ".RData", sep = "")  
save.image(fileName)
# ------------------------------------------------------------------------------