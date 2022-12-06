load("F5_output.RData")

# load data
load("F4_res_pollenData_repeat_10_times_f2000.RData")

true = Label_[,1]
true_num =  Label_[,2] 

# R packages
library("ggplot2","NMI")
library("ClusterR")
# import helper functions
source('../My_code_myFuns/Add_Function_FULL.R')

res_l20_NMF = opt_out = outs_list[[1]][[8]]
get_NMI_Results1(list(list(res_l20_NMF))) # [1] 0.8678727
# ------------------------------------------------------------------------------
rowSum = apply(opt_out$W, 1, sum)
id = which(rowSum!=0)

W = Wmat = opt_out$W[id,]
H = opt_out$H
geneNames = rownames(W)

# ------------------------------------------------------------------------------
# get w and h
# HH = H[c(2,1,4,3), order(true_num)]
# WW = W[ ,c(2,1,4,3)]

HH = H[, order(true_num)]
WW = W

preLab = apply(HH,2,which.max)
cell_type_mat = data.frame(true=true[order(true_num)],
                           true_num=true_num[order(true_num)],
                           preLab=preLab)

# ------------------------------------------------------------------------------
# z-scores method for module members
zmodules = list() 
for(i in 1:4){
  s = WW[,i]
  zs = (s-mean(s))/sd(s)
  print(length(which(zs>1.5)))
  zmodules[[i]] = list(genes = geneNames[which(abs(zs)>1.5)],
                       cells = cell_type_mat$true[cell_type_mat$preLab==i])
}

# ------------------------------------------------------------------------------
# plot heatmap figure
library('ComplexHeatmap')
library("circlize")

get_class_order = function(WW){
  # order for the row
  # Step 1
  classMember = score = rep(0,nrow(WW)) 
  for(i in 1:nrow(WW)){
    j = which.max(WW[i,])
    classMember[i] = j
    score[i] = WW[i,j]
  }
  # Step 2
  id1 = order(classMember)
  dat = dat_raw = data.frame(gene_id = id1,
                             gene_class = classMember[id1],
                             score =score[id1]) 
  for(i in 1:ncol(WW)){
    rows = which(dat$gene_class==i)
    order_ids = order(-dat$score[rows])
    dat[rows,] = dat[rows[order_ids],]
  }
  return(dat)
}

res_W_order = get_class_order(WW)
res_H_order = get_class_order(t(HH))

col_fun = colorRamp2(c(0, 2), c("white", "red"))
col_fun(seq(0, 2))
ht1 = Heatmap(WW[res_W_order$gene_id,], name=paste("W"), 
              cluster_rows=F, cluster_columns=F, col = col_fun,
              show_row_names=F,show_column_names=F)

col_fun = colorRamp2(c(0, 4), c("white", "red"))
col_fun(seq(0, 4))
ht2 = Heatmap(HH[,res_H_order$gene_id], name=paste("H"), 
              cluster_rows=F, cluster_columns=F, col = col_fun,
              show_row_names=F,show_column_names=F)


png(file="F6_W.png",  width = 2.6, height = 7.5, units="in",res= 1200)
print(ht1)
dev.off()

png(file="F6_H.png",  width = 4, height = 2.2, units="in",res= 1200)
print(ht2)
dev.off()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

save(opt_out, pollenData, zmodules, cell_type_mat, file = "F6_outputData.RData")