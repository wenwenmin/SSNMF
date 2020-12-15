# ------------------------------------------------------------------------------
# remove all variables
rm(list = ls())
# ------------------------------------------------------------------------------
# R packages
library("ggplot2","NMI","sparcl")
library("NMI")
library("sparcl")
library("ClusterR") # for external_validation function
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
# ------------------------------------------------------------------------------
# parameter setting
k = 120 
rho = .1 
gama = 4
palm_tol = 1e-3 
rho_tol = 1e-3
palm_maxiter = 200
rho_maxiter = 12

# import helper functions
source('Function/Add_Function_FULL.R')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# test different methods
repeatNum = 10
rawResMat = matrix(0, 12, repeatNum)  
for(j in 1:repeatNum){
  # initialize 
  print(j)
  W0H0 = get_unif_init(X, ncl, nstart0=j)
  Winit = W0H0$W; Hinit = W0H0$H
  
  source('Function/NMF_BCD.R')
  out0 = NMF(X, Winit, Hinit, tol=1e-3, maxit=500)
  NMI0 = get_NMI_Results(list(list(out0)))
  
  source('Function/ONMF_new.R')
  out1 = ONMF(X, Winit, Hinit, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
  NMI1 = get_NMI_Results(list(out1))
  
  source('Function/SONMF_L20_new.R')
  out2 = SONMF_L20(X, Winit, Hinit, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
  NMI2 = get_NMI_Results(list(out2))
  
  source('Function/SONMF_L20_FixRho_new.R')
  ou3 = SONMF_L20_FixRho(X, Winit, Hinit, k, rho=1, palm_tol=1e-3, palm_maxiter)
  NMI3 = get_NMI_Results(list(list(ou3)))
  
  source('Function/SONMF_CL0_new.R')
  ou4 = SONMF_CL0(X, out0$W, out0$H, k, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
  NMI4 = get_NMI_Results(list(ou4))
  
  source('Function/SONMF_L0_new.R')
  ou5 = SONMF_L0(X, out0$W, out0$H, k*3, rho, gama, palm_tol, rho_tol, palm_maxiter, rho_maxiter)
  NMI5 = get_NMI_Results(list(ou5))
  
  source('Function/SNMF_L0_new.R')
  out6 = SNMF_L0(X, out0$W, out0$H, k*3, palm_tol, palm_maxiter)
  NMI6 = get_NMI_Results(list(list(out6)))
  
  source('Function/SNMF_L20_new.R')
  out7 = SNMF_L20(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
  NMI7 = get_NMI_Results(list(list(out7)))
  
  source('Function/SNMF_CL0_new.R')
  out8 = SNMF_CL0(X, out0$W, out0$H, k, palm_tol, palm_maxiter)
  NMI8 = get_NMI_Results(list(list(out8)))
  
  # K-means
  # set.seed(j)
  out9 = kmeans(t(X), centers=t(Winit), nstart = 1)
  # out9 = kmeans(t(X), t(Winit))
  NMI9 = get_NMI(true,out9$cluster)
  
  # sparse K-means with L0 penalty
  source('Function/sparse_kmeans.R')
  set.seed(j)
  out10 = L0_kmeans(t(X), ncl, k, nstart = 1)
  NMI10 = get_NMI(true,out10$cluster)
  
  # sparse K-means with L1 penalty
  # library("sparcl")
  # set.seed(j)
  r11 = KMeansSparseCluster(t(X), centers=t(Winit), wbounds=5.6, nstart = 10)
  NMI11 = get_NMI(true,r11[[1]]$Cs)
  NMI11
  
  rawResMat[,j] = c(NMI0, NMI1, NMI2, NMI3, NMI4, NMI5, NMI6, NMI7, NMI8, NMI9, NMI10, NMI11)
  
  # display result for each testing
  NMIs_diffMet = rawResMat[,j]
  names(NMIs_diffMet) = c("NMF", "ONMF","ONMF_L20", "ONMF_L20_rho","ONMF_CL0","ONMF_L0",
                          "NMF_L0","NMF_L20","NMF_CL0","Kmeans","Kmeans_L0","Kmeans_L1")
  print(as.matrix(NMIs_diffMet))
}
row.names(rawResMat) = c("NMF", "ONMF","ONMF_L20", "ONMF_L20_rho","ONMF_CL0","ONMF_L0",
                         "NMF_L0","NMF_L20","NMF_CL0","Kmeans","Kmeans_L0","Kmeans_L1")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# display results
ResMat = rawResMat # for backup

# sorting
ResMat = ResMat[order(-apply(ResMat,1,mean)),]

# NMI
df_names = row.names(ResMat) 
df = data.frame(method=factor(x=df_names, levels=df_names), 
                avg=apply(ResMat,1,mean), sd =apply(ResMat,1,sd))
print(df)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# plot figure
fig2 = ggplot(df, aes(x=method, y=avg)) + #theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  #geom_bar(stat="identity", color="black", fill="gray", size=0.62, width=0.6) +
  geom_bar(stat="identity", aes(color=method, fill=method), size=0.62, width=0.6) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2, size=1) +
  theme(axis.text = element_text(colour = "black", size = rel(1.5))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"))+
  theme(axis.line  = element_line(size = rel(1.5), colour = "black")) +
  theme(axis.ticks = element_line(size = rel(1.5), colour = "black")) +
  theme(axis.title.y = element_text(size = rel(1.5))) +
  theme(legend.position = "none") + labs(title=NULL, x =NULL, y = "NMI")

fig2 + scale_x_discrete(labels=c( "NMF" = "NMF", 
                           "ONMF" = "ONMF",
                           "ONMF_L20" = "ONMF_L20", 
                           "ONMF_L20_rho" = expression(ONMF_L20_~rho),
                           "ONMF_CL0" = "ONMF_CL0",
                           "ONMF_L0" = "ONMF_L0",
                           "NMF_L0" = "NMF_L0",
                           "NMF_L20" = "NMF_L20",
                           "NMF_CL0" = "NMF_CL0",
                           "Kmeans" = "Kmeans",
                           "Kmeans_L0" = "Kmeans_L0",
                           "Kmeans_L1" = "Kmeans_L1"))

ggsave("Figure/res_barplot.pdf", width = 8, height = 8, units = "in")
ggsave("Figure/res_barplot.png", width = 8, height = 8, units = "in")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# get time
time = proc.time() - ptm; print(time)
fileName = paste("Result/Result_simulation", ".RData", sep = "")  
save.image(fileName)
################################################################################