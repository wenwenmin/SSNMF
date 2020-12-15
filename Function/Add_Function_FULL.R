##############################################################
# Add functions
##############################################################
# get bar plot figure
get_barplot = function(ResMat, type="NMI"){
  # ----------------------------------
  # ResMat is a matrix
  # its each row denotes a method
  # its each column denotes a repeat
  # ----------------------------------
  df1 = data.frame(avg=apply(ResMat,1,mean), sd =apply(ResMat,1,sd))
  # sorting
  id = order(df1$avg, decreasing=T)
  df = data.frame(method=factor(x=row.names(df1)[id], levels=row.names(df1)[id]), 
                  avg = df1$avg[id], sd = df1$sd[id])
  fig = ggplot(df, aes(x=method, y=avg)) + #theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    #geom_bar(stat="identity", color="black", fill="gray", size=0.62, width=0.6) +
    geom_bar(stat="identity", aes(color=method, fill=method), size=0.62, width=0.6) +
    geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2, size=1) +
    theme(axis.text = element_text(colour = "black", size = rel(1.5))) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"))+
    # X和Y轴的线设
    theme(axis.line  = element_line(size = rel(1.5), colour = "black")) +
    # 坐标轴上的小线设
    theme(axis.ticks = element_line(size = rel(1.5), colour = "black")) +
    theme(axis.title.y = element_text(size = rel(1.5))) +
    theme(legend.position = "none") + labs(title=NULL, x =NULL, y = type)
  
  ggsave(paste("F2_", type, ".pdf", sep = ""), width = 7, height = 7, units = "in")
  return(list(df=df1, df_order=df, fig=fig))
}

# true_cluster = get_true_num_label(true)
get_true_num_label = function(true){
  clusterNames = unique(true)
  true_num_label = rep(0,length(true))
  i = 1
  for(c in clusterNames){
    true_num_label[true==c] = i
    i = i + 1
  }
  return(true_num_label)
}

# Initialization method 1
get_unif_init = function(X, ncl, nstart0=100){
  d = nrow(X)
  n = ncol(X)
  set.seed(nstart0)
  W = matrix(runif(d*ncl, 0,1), d, ncl)
  set.seed(nstart0*2)
  H = matrix(runif(ncl*n, 0,1), ncl, n)
  return(list(W=W,H=H))
}

# Initialization method 2
get_kmeans_init = function(X, ncl, nstart0=100){
  # X : row is feature and column is sample
  print(nstart0)
  r0 = kmeans(t(X), ncl, nstart = nstart0)
  W = t(r0$centers)
  H = matrix(0,ncl,dim(X)[2])
  labs = r0$cluster
  for(i in 1:ncl){
    H[i,which(labs==i)] = 1
  }
  return(list(W=W,H=H,cluster=r0$cluster))
}
# (1) Eg. from NetSONMF output import opt W and H
get_optW_optH = function(out,i=NULL){
  if(is.null(i)){
    optW = out[[length(out)]]$W
    optH = out[[length(out)]]$H
    colnames(optH) = 1:ncol(optH)}else{
      optW = out[[i]]$W
      optH = out[[i]]$H
      colnames(optH) = 1:ncol(optH)
    }
  # Prediction label
  preLab = apply(optH,2,which.max); names(preLab) = NULL
  return(list(preLab=preLab, optW=optW, optH=optH))
}
##############################################################
# (2) get NMI score
get_NMI = function(true,preLab){
  t1 = data.frame(1:length(true), true)
  t2 = data.frame(1:length(true), preLab)
  return(NMI_v2(t1,t2)$value)
}
##############################################################
# (3) NMI scores for different methods
get_NMI_Results = function(Results, opt_i=NULL, m="nmi"){
  NMI_Res = c()
  for(i in 1:length(Results)){
    temp = Results[[i]]
    Labs = get_optW_optH(temp, opt_i)$preLab
    NMI_Res[i] = external_validation(true, Labs, method=m)
  }
  names(NMI_Res) =  names(Results) 
  return(NMI_Res)
}

get_NMI_Results1 = function(Results, opt_i=NULL){
  NMI_Res = c()
  for(i in 1:length(Results)){
    temp = Results[[i]]
    Labs = get_optW_optH(temp, opt_i)$preLab
    NMI_Res[i] = get_NMI(true,Labs)
  }
  names(NMI_Res) =  names(Results) 
  return(NMI_Res)
}

get_purity_Results = function(Results, true_lab, opt_i=NULL){
  NMI_Res = c()
  for(i in 1:length(Results)){
    temp = Results[[i]]
    Labs = get_optW_optH(temp, opt_i)$preLab
    NMI_Res[i] = external_validation(true_lab, Labs, method="purity")
  }
  names(NMI_Res) =  names(Results) 
  return(NMI_Res)
}

# print(external_validation(true_cluster, Labs, method = "nmi"))
##############################################################
# (4) Network Enrichment of gene set
NetworkEnrichment.of.module = function(PPINetwork,mgIDs){
  # module gene IDs
  B = PPINetwork
  geneNum = dim(B)[1]
  Num.mgIDs.edges = sum(B[mgIDs,mgIDs])/2
  q = Num.mgIDs.edges-1
  m = sum(B)/2
  n = geneNum*(geneNum-1)/2-m
  k = length(mgIDs)*(length(mgIDs)-1)/2
  
  PValue = phyper(q, m, n, k, lower.tail = F)
  FC = (Num.mgIDs.edges/k)/(m/(m+n)) 
  return(c(length(mgIDs),Num.mgIDs.edges,FC,PValue))
}
##############################################################
# (6) Network Enrichment of gene set for different methods
get_geneEdge_Results = function(Results, opt_i=NULL){
  Res = c()
  for(i in 1:length(Results)){
    temp = Results[[i]]
    optW = get_optW_optH(temp, opt_i)$optW
    mgIDs = which(apply(optW*optW,1,sum)!=0)
    tempRes = NetworkEnrichment.of.module(A,mgIDs)   
    Res = cbind(Res,tempRes)
  }
  colnames(Res) =  names(Results) 
  row.names(Res)=  c("#gene","#edge","FC","p_value")
  return(Res)
}
##############################################################
# The NMI function from "NMI" R package
# This Function Computes the Normalized Mutual Information of Community Strctures in Network
NMI_v2 <- function(X,Y)
{
# X and Y are two matrics representing two partions
# whose first variables are the node id and the second variables are cluster index
  
# Convert X and Y to matrix (mainly for data frames)
X<-as.matrix(X)
Y<-as.matrix(Y)

# Split both X and Y to clusters
XC<-lapply(split(X[,1],X[,2] ), matrix, ncol=1)
YC<-lapply(split(Y[,1],Y[,2] ), matrix, ncol=1)

# Set up the probability matrix
P<-matrix(NA, nrow=length(XC),ncol=length(YC))
for(i in 1:dim(P)[1])
{
  for(j in 1:dim(P)[2])
  {
    P[i,j]<-length(which(as.vector(unlist(XC[i])) %in% as.vector(unlist(YC[j]))))
  }
}
# PI
Pi<-rowSums(P)
# PJ
Pj<-colSums(P)
# P
PP<-sum(P)
# Numerator
D<-0
for(i in 1:dim(P)[1])
{
  for(j in 1:dim(P)[2])
  {
    if(Pi[i]*Pj[j]==0){D<-D+0}
    else{
    temp<-P[i,j]*PP/Pi[i]/Pj[j]
    if(temp==0){D<-D+0}
    if(temp!=0){D<-D+P[i,j]*log(temp)}}
  }
}
D<--2*D
# Denominator 1
N1<-0
for(i in 1:dim(P)[1])
{
  if(PP==0){N1<-N1+0}
  else{
  temp<-Pi[i]/PP
  if(temp==0){N1<-N1+0}
  if(temp!=0){N1<-N1+Pi[i]*log(temp)}}
  
}
# Denominator 2
N2<-0
for(j in 1:dim(P)[2])
{
  if(PP==0){N2<-N2+0}
  else{
  temp<-Pj[j]/PP
  if(temp==0){N2<-N2+0}
  if(temp!=0){N2<-N2+Pj[j]*log(temp)}}
}

if(N1+N2==0){NMI<-0}
if(N1+N2!=0){NMI<-D/(N1+N2)}

# Retuen
object<-list(value=NMI)
object
}
##############################################################