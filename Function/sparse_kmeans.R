################################################################################
################################################################################
L0_kmeans = function(X, ncl, TopK = 10, nstart=10, maxiter=6){
  # ----------------------------------------------------------------------------
  # >> l0km = L0Kmeans(X, ncl, TopK = 2000)
  # X is the data, nxp
  # K is the number of clusters desired
  # Start with equal weights on each feature
  # centers is a kxp matrix
  # ----------------------------------------------------------------------------
  w0 = rep(1, ncol(X)) 
  for(kk in 1:maxiter){
    # update data matrix X
    X_ = X[, which(w0==1)]
    km = kmeans(X_, ncl, nstart)
    cluster = km$cluster
    # print(get_NMI(true,cluster))
	
    # calculate the 'a' vector
    a = calculate_a_vector(X, cluster)
    # update w vector 
    w = 0*w0
    w[order(a,decreasing=T)[1:TopK]] = 1
    if(sum(abs(w-w0))/sum(abs(w0))<1e-4){
      break
    }
    w0 = w
  }
  # # re-run Kmeans with the final w
  # X_ = X[, which(w==1)]
  # km = kmeans(X_, ncl, nstart = nstart0)
  # cluster = km$cluster
  return(list(w=w, cluster=cluster))
}
# ------------------------------------------------------------------------------
calculate_a_vector  = function(X,cluster){
  # X is the data, nxp with samples by features
  n = nrow(X)
  p = ncol(X)
  # initialize a with zero vector
  a = rep(0, p)
  for(j in 1:p){
    X_j = X[,j]
    sum1 = sum(dist(X_j))/n
    
    sum2 = 0
    for(k in unique(cluster)){
      # members in cluster k
      id = which(cluster==k)
      n_k = length(id)
      X_j_k = X_j[id]
      sum2 = sum2 + sum(dist(X_j_k))/n_k
    }
    a[j] = sum1 - sum2
  }
  return(a)
}
################################################################################
################################################################################
L1_kmeans = function(X, ncl, s = 0.1, nstart0=10, maxiter=10){
  # X is the data, nxp with samples by features
  # ncl is the number of clusters
  w0 = rep(1, ncol(X)) 
  for(kk in 1:maxiter){
    # update data matrix X: multiply each column by a constant
    X_ =  t(sqrt(w0[w0!=0])*t(X[,w0!=0]))
    # X_ =  sweep(X[,w0!=0], 2, sqrt(w0[w0!=0]), "*")
    km = kmeans(X_, ncl, nstart = nstart0)
    cluster = km$cluster
    # print(get_NMI(true,cluster))
    
    # calculate the 'a' vector
    a = calculate_a_vector(X, cluster)
    
    # update w vector 
    w = L1L2_Proj(a, s)
    if(sum(abs(w-w0))/sum(abs(w0))<1e-4){
      break
    }
    w0 = w
  }
  return(list(w=w, cluster=cluster))
}
# ------------------------------------------------------------------------------
L1L2_Proj = function(y, c0){
  # min_x ||x-y||_2^2 s.t. ||x||_1<c and ||x||_2=1
  # ICLR 2013 Block Coordinate Descent for Sparse NMF [ncnmf]
  # CJE 2017 A Novel Sparse Penalty for Singular Value Decomposition (Algorithm 4)
  # c0 in [0,1] sparse level
  
  c = (sqrt(length(y))-1)*c0 + 1; # print(c)
  
  v = abs(y)
  if(c < 1){
    print("Full sparsity")
    return(0*y)
  }
  if(c == 1){
    ID = which.max(abs(y))
    y = 0*y
    y[ID] = 1
    return(y)
  }
  if(sum(v)^2/sum(v^2) < c^2){
    print("No sparsity")
    return(y/sqrt(sum(y^2)))
  }
  
  z = v[order(v, decreasing = F)]
  p = length(z)
  for(i in 1:p){
    #print(i)
    eta = z[i]
    zz = z - eta;
    #if(sum(zz[zz>0])^2/sum(zz[zz>0]^2) < c^2){
    if(sum(zz[zz>0])^2 < c^2*sum(zz[zz>0]^2)){
      k_star = p - i + 1
      break
    }
  }
  z_top = z[(p-k_star+1):p]
  s1 = sum(z_top)/k_star
  s2 = sum(z_top^2)/k_star
  
  eta = s1 - sqrt(c^2*(s2-s1^2)/(k_star-c^2)) #4~5
  x = v-eta; x[x<0] = 0
  x = sign(y)*(x)
  x = x/sqrt(sum(x^2))
  return(x)
}  
################################################################################
################################################################################