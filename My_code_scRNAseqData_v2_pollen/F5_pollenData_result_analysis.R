# load data
load("F4_res_pollenData_repeat_10_times_f2000.RData")

true = Label_[,1]
true_num =  Label_[,2] 

# R packages
library("ggplot2","NMI")
library("ClusterR")
# import helper functions
source('../My_code_myFuns/Add_Function_FULL.R')

# load("../My_data/1_pollenData/pollenData.RData")
# dataSet = pollenData
# true = dataSet$sample.group
# true2 = get_true_num_label(dataSet$sample.group)

NMI_res_mat = matrix(0, 11, repeatNum) 
row.names(NMI_res_mat) = c("NMF", "ONMF","ONMF_L20", "ONMF_L20_Rho","ONMF_CL0","ONMF_L0",
                           "NMF_L0","NMF_L20","NMF_CL0","Kmeans","L0_Kmeans")
purity_res_mat = NMI_res_mat

for(j in 1:length(outs_list)){
  out = outs_list[[j]]
  out0 = out[[1]]
  out1 = out[[2]]
  out2 = out[[3]]
  out3 = out[[4]]
  out4 = out[[5]]
  out5 = out[[6]]
  out6 = out[[7]]
  out7 = out[[8]]
  out8 = out[[9]]
  out9 = out[[10]]
  out10 = out[[11]]
  
  NMI0 = get_NMI_Results1(list(list(out0)))
  NMI1 = get_NMI_Results1(list(out1))
  NMI2 = get_NMI_Results1(list(out2))
  NMI3 = get_NMI_Results1(list(list(out3)))
  NMI4 = get_NMI_Results1(list(out4))
  NMI5 = get_NMI_Results1(list(out5))
  NMI6 = get_NMI_Results1(list(list(out6)))
  NMI7 = get_NMI_Results1(list(list(out7)))
  NMI8 = get_NMI_Results1(list(list(out8)))
  NMI9 = get_NMI(true,out9$cluster) 
  NMI10 = get_NMI(true,out10$cluster)
  
  NMI_res_mat[,j] = c(NMI0, NMI1, NMI2, NMI3, NMI4, NMI5, NMI6, NMI7, NMI8, NMI9, NMI10)
  
  
  purity0 = get_purity_Results(list(list(out0)), true_num)
  purity1 = get_purity_Results(list(out1), true_num)
  purity2 = get_purity_Results(list(out2), true_num)
  purity3 = get_purity_Results(list(list(out3)), true_num)
  purity4 = get_purity_Results(list(out4), true_num)
  purity5 = get_purity_Results(list(out5), true_num)
  purity6 = get_purity_Results(list(list(out6)), true_num)
  purity7 = get_purity_Results(list(list(out7)), true_num)
  purity8 = get_purity_Results(list(list(out8)), true_num)
  purity9 = external_validation(true_num, out9$cluster, method="purity")
  purity10= external_validation(true_num, out10$cluster, method="purity") 
  
  purity_res_mat[,j] = c(purity0, purity1, purity2, purity3, purity4, purity5, 
                         purity6, purity7, purity8, purity9, purity10)
}
# ------------------------------------------------------------------------------

################################################################################
# gcR: Get Compare Results with different evaluation metrics
gcr = function(outs_list, eval_met = "purity"){
  # step1: get different methods predicted cluster label 
  pre_label_list = c()
  for(j in 1:length(outs_list)){
    out = outs_list[[j]]
    # each method each element
    pre_list = c()
    pre_list[[1]] = get_optW_optH(list(out[[1]]),i=NULL)$preLab # NMF
    pre_list[[2]] = get_optW_optH(out[[2]],i=NULL)$preLab       # ONMF  
    pre_list[[3]] = get_optW_optH(out[[3]],i=NULL)$preLab       # SONMF_L20_new 
    pre_list[[4]] = get_optW_optH(list(out[[4]]),i=NULL)$preLab # SONMF_L20_FixRho
    pre_list[[5]] = get_optW_optH(out[[5]],i=NULL)$preLab       # SONMF_CL0
    pre_list[[6]] = get_optW_optH(out[[6]],i=NULL)$preLab       # SONMF_L0
    pre_list[[7]] = get_optW_optH(list(out[[7]]),i=NULL)$preLab # SNMF_L0
    pre_list[[8]] = get_optW_optH(list(out[[8]]),i=NULL)$preLab # SNMF_L20
    pre_list[[9]] = get_optW_optH(list(out[[9]]),i=NULL)$preLab # SNMF_CL0
    pre_list[[10]] = out[[10]]$cluster # K-means
    pre_list[[11]] = out[[11]]$cluster # L0 K-means
    pre_label_list[[j]] = pre_list
  }
  # step2:
  mat = matrix(0, 11, 10) 
  for(j in 1:length(pre_label_list)){
    pre_list = pre_label_list[[j]]
    for(i in 1:length(pre_list))
      mat[i,j] = external_validation(true_num, pre_list[[i]], method = eval_met)
  }
  row.names(mat) = c("NMF", "ONMF","ONMF_L20", "ONMF_L20_Rho","ONMF_CL0","ONMF_L0",
                     "NMF_L0","NMF_L20","NMF_CL0","Kmeans","L0_Kmeans")
  
  avg_sd = data.frame(avg=apply(mat,1,mean), sd=apply(mat,1,sd))
  return(avg_sd)
}
################################################################################
# display results
NMI_avg_sd = data.frame(avg=apply(NMI_res_mat,1,mean), sd=apply(NMI_res_mat,1,sd))
purity_avg_sd = data.frame(avg=apply(purity_res_mat,1,mean), sd=apply(purity_res_mat,1,sd))
purity_avg_sd2 = gcr(outs_list, eval_met = "purity")
entropy_avg_sd = gcr(outs_list, eval_met = "entropy")

# ------------------------------------------------------------------------------
# plot bar plot
res_plot_NMI = get_barplot(NMI_res_mat,"NMI")
res_plot_purity = get_barplot(purity_res_mat,"purity")

# ------------------------------------------------------------------------------
# gct: Get a Column of Table
gct = function(NMI_avg_sd, na="NMI"){
  NMI_sd = c()
  for(i in 1:nrow(NMI_avg_sd)){
    
    t1 = format(round(NMI_avg_sd[i,1]*100,2), nsmall = 2)
    t2 = format(round(NMI_avg_sd[i,2]*100,2), nsmall = 2)
    
    NMI_sd[i] = paste(t1,"+",t2)
  }
  dat = data.frame(row.names = row.names(NMI_avg_sd), NMI_sd=NMI_sd)
  names(dat) = na
  return(dat)
}

# get Table 1
# ------------------------------------------------------------------------------
# The second column
c0 = data.frame(row.names=row.names(NMI_res_mat), 
                gene=c("all","all","2000","2000","XX","XX","XX","2000","XX","all","2000"))
names(c0) = "gene"

c1 = gct(NMI_avg_sd, "NMI+sd")
c2 = gct(purity_avg_sd, "purity+sd")
c3 = gct(entropy_avg_sd, "entropy+sd")
TabData =cbind(cbind(c0,cbind(c1,c2)),c3)

library(xtable)
pollenDataRes = xtable(TabData)
save(pollenDataRes, file = "F5_pollenDataRes.RData")

# ------------------------------------------------------------------------------
#save(NMI_res_mat, purity_res_mat, file = "F5_output.RData")