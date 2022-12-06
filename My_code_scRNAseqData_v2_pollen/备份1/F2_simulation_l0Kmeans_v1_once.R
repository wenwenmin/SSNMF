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
# simulation function
sim = function(mu, f=500, u=0.9, nstr=2){
  set.seed(nstr)
  D = u*matrix(rnorm(60*f), f ,60)
  D[1:60, 1:20] = matrix(rnorm(60*20, mean = mu, sd = 1), 60 ,20)
  D[31:90, 21:40] = matrix(rnorm(60*20, mean = mu, sd = 1), 60 ,20)
  D[61:120, 41:60] = matrix(rnorm(60*20, mean = mu, sd = 1), 60 ,20)
  # Non-negative condition
  D = abs(D)
  return(D)
}
# generate a dataset
X = sim(mu=1, f=500, u=0.9, nstr=2)
# vector of true cluster labels
ncl = 3; true = rep(1:3, each=20)  
# ------------------------------------------------------------------------------
# 导入辅助函数
source('../My_code_myFuns/Add_Function_FULL.R')
# ------------------------------------------------------------------------------
# 初始化条件
# W0H0 = get_kmeans_init(X,ncl,nstart0=3)
W0H0 = get_unif_init(X, ncl, nstart0=1)
Winit = W0H0$W; Hinit = W0H0$H

# ------------------------------------------------------------------------------
# 参数设置
k = 120 
rho = .1 
gama = 4
palm_tol = 1e-3 
rho_tol = 1e-3
palm_maxiter = 200
rho_maxiter = 12


W0H0 = get_unif_init(X, ncl, nstart0=2)
Winit = W0H0$W; Hinit = W0H0$H

# ------------------------------------------------------------------------------
source('../My_code_myFuns/NMF_BCD.R')
out0 = NMF(X, Winit, Hinit, tol=1e-3, maxit=500)
NMI0 = get_NMI_Results(list(list(out0)))

# ------------------------------------------------------------------------------
source('../My_code_myFuns/ONMF_new.R')
out1 = ONMF(X, Winit, Hinit, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
NMI1 = get_NMI_Results(list(out1))

# ------------------------------------------------------------------------------
source('../My_code_myFuns/SONMF_L20_new.R')
out2 = SONMF_L20(X, Winit, Hinit, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
NMI2 = get_NMI_Results(list(out2))

# ------------------------------------------------------------------------------
source('../My_code_myFuns/SONMF_L20_FixRho_new.R')
ou3 = SONMF_L20_FixRho(X, Winit, Hinit, k, rho=1, palm_tol=1e-3, palm_maxiter)
NMI3 = get_NMI_Results(list(list(ou3)))

# ------------------------------------------------------------------------------
source('../My_code_myFuns/SONMF_CL0_new.R')
ou4 = SONMF_CL0(X, out0$W, out0$H, 60, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
NMI4 = get_NMI_Results(list(ou4))

# ------------------------------------------------------------------------------
source('../My_code_myFuns/SONMF_L0_new.R')
ou5 = SONMF_L0(X, out0$W, out0$H, 120*3, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
NMI5 = get_NMI_Results(list(ou5))

# ------------------------------------------------------------------------------
source('../My_code_myFuns/SNMF_L0_new.R')
out6 = SNMF_L0(X, out0$W, out0$H, k*3, palm_tol, palm_maxiter)
NMI6 = get_NMI_Results(list(list(out6)))

# ------------------------------------------------------------------------------
source('../My_code_myFuns/SNMF_L20_new.R')
out7 = SNMF_L20(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
NMI7 = get_NMI_Results(list(list(out7)))

# ------------------------------------------------------------------------------
source('../My_code_myFuns/SNMF_CL0_new.R')
out8 = SNMF_CL0(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
NMI8 = get_NMI_Results(list(list(out8)))

# ------------------------------------------------------------------------------
# out9 = kmeans(t(X), ncl, nstart = 1)
# centers with initialization
out9 = kmeans(t(X), t(Winit))
NMI9 = get_NMI(true,out9$cluster)

source('../My_code_myFuns/L0_kmeans.R')
out10 = L0_kmeans(t(X), ncl, k, nstart0=1)
NMI10 = get_NMI(true,out10$cluster)
# ------------------------------------------------------------------------------

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

print(NMIs_diffMet)
# ------------------------------------------------------------------------------