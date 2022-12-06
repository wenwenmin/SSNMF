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

ncl = length(unique(true))

A = as.matrix(dataSet$ppi_Net)

true = get_true_num_label(dataSet$sample.group)
true_char = dataSet$sample.group
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
for(k in ks){
  ResMat = matrix(0, 3, 10)  
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
    preLab  = apply(out7$H,2,which.max)
    NMI     = get_NMI(true_char, preLab)
    # NMI     = external_validation(true, preLab, method= "nmi")
    purity  = external_validation(true, preLab, method= "purity")
    entropy = external_validation(true, preLab, method= "entropy")
    
    ResMat[,j] = c(NMI, purity, entropy)
    row.names(ResMat) = c("nmi", "purity", "entropy")
  }
  res_out[[t]] = ResMat; t = t + 1

  print(c(k,apply(ResMat,1,mean)))
}

NMIs = c()
for(i in 1:length(res_out)){
  NMIs[i] = mean(res_out[[i]]["nmi",])
}
plot(NMIs)

time = proc.time() - ptm; print(time)
save(res_out, NMIs, ks, file = "F11_pollen_NMF_l20.RData")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
library(ggplot2)
dat = data.frame(para_c=ks, no_zero_num = NMIs)
mytheme = theme(axis.line=element_line(size=rel(1.2),colour="white"),
                axis.ticks=element_line(size=rel(1.2)),
                plot.title = element_text(hjust = 0.5),
                axis.title = element_text(size = rel(1.3)),
                axis.text = element_text(colour = "black",size=rel(1.1)),
                legend.title =element_text(size=1, colour="white"),
                legend.text=element_text(size=rel(1.3), colour="black"))

p1 = ggplot(dat, aes(x=para_c, y=no_zero_num)) + geom_line() + geom_point(shape = 1, size = 1.5) +
  labs(title = NULL, x = "Number of genes", y = "NMI") + mytheme
p1

ggsave("F11_pollen_NMF_l20.png", width = 4, height = 3, units = "in")
# ------------------------------------------------------------------------------
