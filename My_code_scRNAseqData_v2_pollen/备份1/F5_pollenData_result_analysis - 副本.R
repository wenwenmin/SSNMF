# load data
load("F4_res_pollenData_repeat_10_times_f2000.RData")

# R packages
library("ggplot2","NMI")
library("ClusterR")
# import helper functions
source('../My_code_myFuns/Add_Function_FULL.R')

load("../My_data/1_pollenData/pollenData.RData")
dataSet = pollenData
true = dataSet$sample.group
true2 = get_true_num_label(dataSet$sample.group)

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
  
  
  purity0 = get_purity_Results(list(list(out0)), true2)
  purity1 = get_purity_Results(list(out1), true2)
  purity2 = get_purity_Results(list(out2), true2)
  purity3 = get_purity_Results(list(list(out3)), true2)
  purity4 = get_purity_Results(list(out4), true2)
  purity5 = get_purity_Results(list(out5), true2)
  purity6 = get_purity_Results(list(list(out6)), true2)
  purity7 = get_purity_Results(list(list(out7)), true2)
  purity8 = get_purity_Results(list(list(out8)), true2)
  purity9 = external_validation(true2, out9$cluster, method="purity")
  purity10= external_validation(true2, out10$cluster, method="purity") 
                                
  purity_res_mat[,j] = c(purity0, purity1, purity2, purity3, purity4, purity5, 
                         purity6, purity7, purity8, purity9, purity10)
}

# ------------------------------------------------------------------------------
# display results
# ------------------------------------------------------------------------------
ResMat = purity_res_mat
ResMat = ResMat[order(-apply(ResMat,1,mean)),]

# NMI平均值和方差
df_names = row.names(ResMat) 
df = data.frame(method=factor(x=df_names, levels=df_names), 
                avg=apply(ResMat,1,mean), sd =apply(ResMat,1,sd))
print(df)

dat <- stack(as.data.frame(t(ResMat)))
dat$ind = factor(dat$ind, levels=row.names(ResMat))
fig2 = ggplot(dat, aes(x=ind, y=values, color=ind)) +
  geom_boxplot(size=1.2) +
  theme(axis.line=element_line(size=rel(1.2)),
        axis.ticks=element_line(size=rel(1.2)),
        axis.title = element_text(size = rel(1.3)),
        axis.text = element_text(colour = "black",size=rel(1.2)),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none",
        legend.title =element_text(size=1, colour="white"),
        legend.text=element_text(size=rel(1.3), colour="black")) +
  labs(title=NULL, x =NULL, y = "NMI")
ggsave("F5_2.pdf", width = 7, height = 7, units = "in")