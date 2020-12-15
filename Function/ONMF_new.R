ONMF = function(X, Winit, Hinit, 
                rho = 1e-8, 
                gama = 1.5, 
                palm_tol = 3e-3, 
                rho_tol = 3e-3, 
                palm_maxiter = 500, 
                rho_maxiter = 100){
  #-----------------------------------------------------------------------------
  # Input:
  # X is a d*n matrix with d rows/features and n columns/samples
  #-----------------------------------------------------------------------------
  
  out <- list() 
  for(iter_rho in 1:rho_maxiter){
    out[[iter_rho]] <- ONMF_subproblem(X, Winit, Hinit, rho, palm_tol, palm_maxiter)
    
    # Orthogonality
    W = out[[iter_rho]]$W
    H = out[[iter_rho]]$H 
    H_orth <- diag(apply(H*H,1,sum)^(-1/2))%*%H
    e_orth <- norm(H_orth%*%t(H_orth)-diag(nrow(Hinit)),'F')/nrow(Hinit)^2
    
    # Convergence
    rho_delta <- norm(W-Winit,'F')/norm(Winit,'F') + norm(H-Hinit,'F')/norm(Hinit,'F')
    if(e_orth<=rho_tol && rho_delta<=rho_tol) break

    # update W0, H0 and rho
    Winit <- out[[iter_rho]]$W
    Hinit <- out[[iter_rho]]$H
    rho <- rho*gama
  }
  return(out)
}

ONMF_subproblem = function(X, W0, H0, rho, palm_tol, palm_maxiter){
  # ----------------------------------------------------------------------------
  # X is a d*n matrix with d rows/features and n columns/samples
  # Both u/2||W||_F^2 and v/2||H||_F^2 make the model more robust
  # ----------------------------------------------------------------------------
  
  u <- 0
  v <- 1e-4
  
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