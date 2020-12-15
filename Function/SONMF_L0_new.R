# Sparse ONMF for W with L0 constrant
SONMF_L0 = function(X, Winit, Hinit, ks, 
                    rho = 1e-8, 
                    gama = 1.5, 
                    palm_tol = 3e-3,
                    rho_tol = 3e-3, 
                    palm_maxiter = 500, 
                    rho_maxiter = 100){
  #-----------------------------------------------------------------------------
  # Input:
  # X is a d*n matrix with d rows/features and n columns/samples
  # ks: ||W||_0 <=ks 
  #-----------------------------------------------------------------------------
  out <- list() 
  for(iter_rho in 1:rho_maxiter){
    out[[iter_rho]] <- SONMF_L0_FixRho(X, Winit, Hinit, ks, rho, palm_tol, palm_maxiter)
    
    # 1) Orthogonality
    W = out[[iter_rho]]$W
    H = out[[iter_rho]]$H 
    
    # 2) Residuals of W and H
    H_orth <- diag(apply(H*H,1,sum)^(-1/2))%*%H
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

SONMF_L0_FixRho = function(X, W0, H0, ks, rho, tol=1e-3, maxiter=200){
  
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
    # extract ks elements for W matrix
    vec_W = c(W)
    ids <- order(vec_W, decreasing = T)[1:ks]
    vec_W[-ids] <- 0 
    W = matrix(vec_W, ncol = ncol(W), byrow=FALSE)
    # W = ncol(W)*W/norm(W,"F")
    
    # --------------------------------------------------------------------------
    # update H
    Hessian <- t(W)%*%W + rho*matrix(1, ncol(W), ncol(W));
    # get the Lipschitz constant with spectral norm
    Lh <- norm(Hessian + (v-rho)*diag(nrow(Hessian)),'2') 
    gradH <- -t(W)%*%X + Hessian%*%H0 + (v - rho)*H0
    H <- H0 - (1/Lh)*gradH;
    H[H<0] <- 0
    
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