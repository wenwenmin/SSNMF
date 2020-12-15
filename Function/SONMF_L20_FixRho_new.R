SONMF_L20_FixRho = function(X, W0, H0, k, 
                            rho = 1e-8, 
                            palm_tol=3e-3, 
                            palm_maxiter=500, 
                            initial="NMF"){
  #-----------------------------------------------------------------------------
  # Input:
  # X is a d*n matrix with d rows/features and n columns/samples
  # k: extract top k features
  #-----------------------------------------------------------------------------
  u <- 0
  v <- 1e-4
  #-----------------------------------------------------------------------------
  #print("Initial with NMF") 
  if(initial=="NMF"){
    out = NMF(X, W0, H0)
    W0 = out$W; H0 = out$H
	}
  #-----------------------------------------------------------------------------
  objs <- c()
  for(iter in 1:palm_maxiter){
    # --------------------------------------------------------------------------
    # update W
    HHT <- H0%*%t(H0);
    # Lipschitz constant with spectral norm
    Lw <- norm(HHT + u*diag(nrow(HHT)),'2') 
    gradW <- -X%*%t(H0) + W0%*%HHT + u*W0
    W <- W0-(1/Lw)*gradW
    W[W<0] <- 0
    # extract features: point multiplication, sum by row, sort
    ids <- order(apply(W*W,1,sum), decreasing = T)[1:k] 
    W[-ids, ] <- 0 # Row sparse
    
    # --------------------------------------------------------------------------
    # update H
    Hessian <- t(W)%*%W + rho*matrix(1, ncol(W), ncol(W));
    # Lipschitz constant with spectral norm
    Lh <- norm(Hessian + (v-rho)*diag(nrow(Hessian)),'2') 
    gradH <- -t(W)%*%X + Hessian%*%H0 + (v - rho)*H0
    H <- H0 - (1/Lh)*gradH;
    H[H<0] <- 0
    
    # --------------------------------------------------------------------------
    # calculate the objective function value
    objs[iter] <- norm(X-W%*%H,'F')^2/2 + 
                  rho*sum(apply(H,2,sum)^2)/2 + 
                  u*norm(W,'F')^2/2 + (v-rho)*norm(H,'F')^2/2
    
    # --------------------------------------------------------------------------  
    # stopping condition
    delta <- norm(W-W0,'F')/norm(W0,'F') + norm(H-H0,'F')/norm(H0,'F');
    if(iter>5 && delta <= palm_tol) break 
    
    # --------------------------------------------------------------------------
    # remember previous iteration result
    W0 <- W; H0 <- H;
  }
  return(list(W=W,H=H,rho=rho,objs=objs))
}