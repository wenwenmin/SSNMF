# please see the latest version by opening 'sparse_kmeans.R'
L0_kmeans = function(X, ncl, TopK = 10, nstart=10, maxiter=6){
  # ----------------------------------------------------------------------------
  # >> l0km = L0Kmeans(X, ncl, TopK = 2000)
  # X is the data, nxp
  # K is the number of clusters desired
  # Start with equal weights on each feature
  
  # Example
  # l0km = L0Kmeans(X, ncl, TopK = 2000, nstart=1)
  # print(get_NMI(true,l0km$cluster)) 
  # X_ = X[, which(l0km$ws==1)]
  # km = kmeans(X_, ncl, nstart = 1)
  # print(get_NMI(true,km$cluster)) 
  # ----------------------------------------------------------------------------
  ws = rep(1, ncol(X)) 
  ws.old <- rnorm(ncol(X))
  niter = 1
  while((sum(abs(ws-ws.old))/sum(abs(ws.old)))>1e-4 && niter<maxiter){
    X_ = X[, which(ws==1)]
    km = kmeans(X_, ncl, nstart)
    cluster = km$cluster
    a = get_featureImportance_score(X,cluster)
    ws = 0*ws
    ws[order(a,decreasing=T)[1:TopK]] = 1
    niter = niter + 1
  }
  out = list(ws=ws, cluster=cluster)
  return(out)
}
# ------------------------------------------------------------------------------
get_featureImportance_score = function(X,cluster){
  # X is the data, nxp with n samples by p features
  a = c()
  for(j in 1:ncol(X)){
    X_j = X[,j]
    sum1 = sum(dist(X_j))
    
    sum2 = 0
    for(c in unique(cluster)){
      X_j_c = X_j[which(cluster==c)]
      sum2 = sum2 + sum(dist(X_j_c))
    }
    a_j = sum1 - sum2
    a = c(a,a_j)
  }
  return(a)
}
# ------------------------------------------------------------------------------