# Sparse ONMF for W with column L0 constrant
SONMF_CL0 = function(X, Winit, Hinit, 
                     k, 
                     rho = 1e-8, 
                     gama = 1.5, 
                     palm_tol = 3e-3,
                     rho_tol = 3e-3, 
                     palm_maxiter = 500, 
                     rho_maxiter = 100){
  #-----------------------------------------------------------------------------
  # Input:
  # X is a d*n matrix with d rows/features and n columns/samples
  # k: ||W_i||_0 <=k for all i 
  #-----------------------------------------------------------------------------
  
  out <- list() 
  for(iter_rho in 1:rho_maxiter){
    out[[iter_rho]] <- SONMF_CL0_FixRho(X, Winit, Hinit, k, rho, palm_tol, palm_maxiter)
   
    # 1) Orthogonality
    W = out[[iter_rho]]$W
    H = out[[iter_rho]]$H 
    
    # 2) Residuals of W and H
    H_orth <- diag((apply(H*H,1,sum)+0.0001)^(-1/2))%*%H
    e_orth <- norm(H_orth%*%t(H_orth)-diag(nrow(H_orth)),'F')/nrow(H_orth)^2
    rho_delta <- norm(W-Winit,'F')/norm(Winit,'F') + norm(H-Hinit,'F')/norm(Hinit,'F')
    
    # 3) cat("SONMF: stopping update rho.\n")
    if(e_orth<=rho_tol && rho_delta<=rho_tol) break
    
    # 4) update parameters
    Winit <- out[[iter_rho]]$W
    Hinit <- out[[iter_rho]]$H
    rho   <- rho*gama
  }
  return(out)
}

SONMF_CL0_FixRho = function(X, W0, H0, k, rho, tol=1e-3, maxiter=200){
  u <- 0
  v <- 1e-4
  objs <- c()  
  for(iter in 1:maxiter){
    # --------------------------------------------------------------------------
    # update W
    HHT <- H0%*%t(H0);
    # Lipschitz constant
    Lw <- norm(HHT + u*diag(nrow(HHT)),'2') 
    gradW <- -X%*%t(H0) + W0%*%HHT + u*W0
    W <- W0-(1/Lw)*gradW
    W[W<0] <- 0
    # extract k features for W[,j]
    for(j in 1:ncol(W)){
      w <- W[,j]
      sort_index <- order(w, decreasing = T)
      ids <- sort_index[1:k]
      w[-ids] <- 0 
      # w = w/sqrt(sum(w^2))  #good?
      W[,j] <- w
    }
    W = ncol(W)*W/norm(W,"F")
    
    # --------------------------------------------------------------------------
    # update H
    Hessian <- t(W)%*%W + rho*matrix(1, ncol(W), ncol(W));
    # get the Lipschitz constant with spectral norm
    Lh <- norm(Hessian + (v-rho)*diag(nrow(Hessian)),'2') 
    gradH <- -t(W)%*%X + Hessian%*%H0 + (v - rho)*H0
    H <- H0 - (1/Lh)*gradH;
    H[H<0] <- 0
    # Standardization Unit
    # H <- diag(apply(H*H,1,sum)^(-1/2))%*%H
    # W <- W%*%diag(apply(H*H,1,sum)^(1/2))
   
     # --------------------------------------------------------------------------
    # objective function value
    objs[iter] <- norm(X-W%*%H,'F')^2/2 + u*norm(W,'F')^2/2 + (v-rho)*norm(H,'F')^2/2
    
    # stopping condition
    delta <- norm(W-W0,'F')/norm(W0,'F') + norm(H-H0,'F')/norm(H0,'F');
    if(iter>10 && delta <= tol) break # cat("stopping update W and H\n")
    
    # remember previous iteration result
    W0 <- W; H0 <- H;
  }
  return(list(W=W,H=H,rho=rho,objs=objs))
}