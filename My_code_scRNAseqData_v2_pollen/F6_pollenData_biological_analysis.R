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

################################################################################
################################################################################
# 1 get heatmap
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
################################################################################
get_heatmap = function(W, H){
  
  res_W_order = get_class_order(W)
  res_H_order = get_class_order(t(H))
  
  col_fun = colorRamp2(c(0, 2), c("white", "red"))
  col_fun(seq(0, 2))
  ht1 = Heatmap(W[res_W_order$gene_id,], name=paste("W"), 
                cluster_rows=F, cluster_columns=F, col = col_fun,
                show_row_names=F,show_column_names=F)
  
  col_fun = colorRamp2(c(0, 4), c("white", "red"))
  col_fun(seq(0, 4))
  ht2 = Heatmap(H[,res_H_order$gene_id], name=paste("H"), 
                cluster_rows=F, cluster_columns=F, col = col_fun,
                show_row_names=F,show_column_names=F)
  
  return(list(ht1=ht1,ht2=ht2))
}
################################################################################
################################################################################
library('ComplexHeatmap')
library("circlize")

figs = get_heatmap(W, H)

png(file="F6_W.png",  width = 2.6, height = 7.5, units="in",res= 1200)
print(figs$ht1)
dev.off()

png(file="F6_H.png",  width = 4, height = 2.2, units="in",res= 1200)
print(figs$ht2)
dev.off()

################################################################################
# 2 get modules
getModules = function(W, H, Label_){
  # get sample names
  true = Label_[,1]
  preLab = apply(H, 2, which.max)
  
  cell_type_mat = data.frame(id=1:ncol(H), true=true, rank=preLab)
  cell_type_mat_order = cell_type_mat[order(preLab),]
  
  modules = list() 
  module_size_mat = data.frame(num_gene=rep(NA,ncol(W)),
                               num_cell=rep(NA,ncol(W)))
  combine_genes = c()
  geneNames = rownames(W)
  for(i in 1:ncol(W)){
    # z-scores method for a gene set
    s = W[,i]
    zs = (s-mean(s))/sd(s)
    
    gene_set = geneNames[which(zs>1.5)]
    # combine_genes = c(combine_genes,which(zs>1.5))
    combine_genes = c(combine_genes,gene_set)
    
    cell_set = cell_type_mat$true[cell_type_mat$rank==i]
    
    modules[[i]] = list(genes=gene_set, cells=cell_set)
    
    module_size_mat[i,1] = length(gene_set)
    module_size_mat[i,2] = length(cell_set)
  }
  return(list(modules = modules, 
              combine_genes = combine_genes,
              module_size_mat = module_size_mat, 
              cell_type_mat = cell_type_mat,
              cell_type_mat_order = cell_type_mat_order))
}

modres = getModules(W, H, Label_)
length(unique(modres$combine_genes)) # 329

write.table(modres$cell_type_mat_order, file = "Result/F6_cell_type_mat_order.csv", append = FALSE, quote = F, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = "double")

biclusters = modres$modules
write.table(biclusters[[1]]$genes, file = "Result/F6_m1_geneList.txt", append = FALSE, quote = F, sep = "",
            eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F, qmethod = "double")


save(opt_out, pollenData, modres, file = "F6_outputData.RData")