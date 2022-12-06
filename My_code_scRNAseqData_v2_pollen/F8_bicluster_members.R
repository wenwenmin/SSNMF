load("F6_outputData.RData")

load("../My_code_KEGG_GOBP_functional_analysis/geneNameData/Result_pollenData_geneNames.RData")

W = opt_out$W
H = opt_out$H

################################################################################
# ------------------------------------------------------------------------------
get_modules_entrezID = function(biclusters, geneNames.symbols.entrezID){
  geneNames_entrezID = geneNames.symbols.entrezID$NCBI.entrezID
  geneNames_symbols = geneNames.symbols.entrezID$symbols
  modules_entrezID = list()
  for(i in 1:length(biclusters)){
    ids = which(geneNames_symbols%in%biclusters[[i]]$genes==T)
    modules_entrezID[[i]] = geneNames_entrezID[ids] 
  }
  return(modules_entrezID)
}
# ------------------------------------------------------------------------------
biclusters = modres$modules
modules = get_modules_entrezID(biclusters, geneNames.symbols.entrezID)

get_cluster_member = function(biclusters,W){
  # Step 1: get total genes
  combine_genes = c()
  for(i in 1:length(biclusters)){
    combine_genes = c(combine_genes, biclusters[[i]]$genes) 
  }
  combine_genes = sort(unique(combine_genes))
  
  # Step 2: get cluster_member
  cluster_member = data.frame(gene = combine_genes, 
                              matrix(0,nrow=length(combine_genes),ncol=length(biclusters))) 
  colnames(cluster_member) = c("gene", paste("factor",1:length(biclusters),sep = "_"))
  
  for(i in 1:length(biclusters)){
    gene_id = biclusters[[i]]$genes
    cluster_member[combine_genes%in%gene_id, i+1] = 1
  }
  # Step 3: summary
  mat = cluster_member
  row.names(mat) = mat$gene
  mat$gene = NULL
 
  W_mat = W[cluster_member$gene,]
  
  cluster_member$level = apply(mat,1,sum)
  return(list(cluster_member=cluster_member,W_mat=W_mat))
}

res = get_cluster_member(biclusters,W)

write.table(res$cluster_member, file = "Result/F8_cluster_member.csv", append = FALSE, quote = F, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, qmethod = "double")

write.table(res$W_mat, file = "Result/F8_cluster_member_weight.csv", append = FALSE, quote = F, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T, qmethod = "double")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Load library
library(gplots)

set_list = list(M1=modules[[1]], 
                M2=modules[[2]], 
                M3=modules[[3]], 
                M4=modules[[4]])
# venn(set_list)
# category_names = c("M1" , "M2 " , "M3", "M4")
# myfilename = 'F9_venn_diagramm.png'
# fig = get_venn_4sets(set_list,myfilename, category_names)
require(VennDiagram)
vp <- venn.diagram(set_list, 
                   fill = 2:5,
                   alpha = 0.3, 
                   lwd = 0.1,
                   lty = 'blank',
                   filename = NULL);
grid.draw(vp)