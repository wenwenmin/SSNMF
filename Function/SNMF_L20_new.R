SNMF_L20 = function(X, W0, H0, k, tol=3e-3, maxiter=500){
  #-----------------------------------------------------------------------------
  # Input:
  # X is a d*n matrix with d rows/features and n columns/samples
  # k: extract top k features
  #-----------------------------------------------------------------------------
  
  # two parameter
  # Both u/2||W||_F^2 and v/2||H||_F^2 make the model more robust
  u <- 0
  v <- 0 #1e-4
  
  objs <- c() # Objective function values
  for(iter in 1:maxiter){
    
    # --------------------------------------------------------------------------
    # update W
    HHT <- H0%*%t(H0);
    # Lipschitz constant with spectral norm
    Lw <- norm(HHT + u*diag(nrow(HHT)),'2') 
    gradW <- -X%*%t(H0) + W0%*%HHT + u*W0
    W <- W0-(1/Lw)*gradW
    W[W<0] <- 0
    # extract features
    ids = order(apply(W*W,1,sum), decreasing = T)[1:k] 
    W[-ids, ] <- 0 # Row sparse
    
    # --------------------------------------------------------------------------
    # update H
    WTW = t(W)%*%W
    # Lipschitz constant
    Lh <- norm(WTW + v*diag(nrow(WTW)),'2') 
    gradH <- -t(W)%*%X + WTW%*%H0 + v*H0
    H <- H0 - (1/Lh)*gradH;
    H[H<0] <- 0

    # --------------------------------------------------------------------------
    # calculate the objective function value
    objs[iter] <- norm(X-W%*%H,'F')^2/2 + u*norm(W,'F')^2/2 + v*norm(H,'F')^2/2
    
    # stopping condition
    delta <- norm(W-W0,'F')/norm(W0,'F') + norm(H-H0,'F')/norm(H0,'F');
    if(iter > 20 && delta <= tol) break
   
    # remember previous iteration result
    W0 <- W; H0 <- H;
  }
  return(list(W=W,H=H,objs=objs))
}