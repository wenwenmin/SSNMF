SONMF_L20_aPALM = function(X, W0, H0, k, 
                           rho = 1e-8, 
                           gama = 1.5, 
                           palm_tol = 3e-3, 
                           rho_tol = 3e-3, 
                           palm_maxiter = 500, 
                           rho_maxiter = 100){
  #-----------------------------------------------------------------------------
  # Input:
  # X is a d*n matrix with d rows/features and n columns/samples
  # k: extract top 'k' features
  #-----------------------------------------------------------------------------
  out <- list() 
  
  #-----------------------------------------------------------------------------
  # solve a subproblem for each rho
  for(iter_rho in 2:rho_maxiter){
    
    out[[iter_rho]] <- SONMF_L20_FixRho_aPALM(X, W0, H0, k, rho, palm_tol, palm_maxiter)
    
    
    # 1) Orthogonality
    W <- out[[iter_rho]]$W
    H <- out[[iter_rho]]$H 
    
    H_orth <- diag(apply(H*H,1,sum)^(-1/2))%*%H
    e_orth <- norm(H_orth%*%t(H_orth)-diag(nrow(H_orth)),'F')/nrow(H_orth)^2
    # 2) Residuals of W and H
    rho_delta <- norm(W-W0,'F')/norm(W0,'F') + norm(H-H0,'F')/norm(H0,'F')
    
    out[[1]]$Orth <- e_orth
    out[[1]]$Residual <- rho_delta
    
    # 3) cat("SONMF: stopping update rho.\n")
    if(e_orth <= rho_tol && rho_delta <= rho_tol) break
    
    # 4) update parameters
    W0 <- out[[iter_rho]]$W
    H0 <- out[[iter_rho]]$H
    rho <- rho*gama
  }
  return(out)
}

SONMF_L20_FixRho_aPALM = function(X, W0, H0, k, rho, palm_tol=1e-3, palm_maxiter=500){
  # ----------------------------------------------------------------------------
  # X is a d*n matrix with d rows/features and n columns/samples
  # Both u/2||W||_F^2 and v/2||H||_F^2 make the model more robust
  # ----------------------------------------------------------------------------
  
  u <- 0
  v <- 1e-4
  
  rw <- 0.5
  Xnrm <- norm(X,'f')
  obj0 <- 0.5*Xnrm*Xnrm
  
  t0 <- 1
  Lw <- 1
  Lh <- 1
  Wm <- W0 
  Hm <- H0 
  
  objs <- c()
  for(iter in 1:palm_maxiter){
    
    # --------------------------------------------------------------------------
    # update W
    Hs <- H0%*%t(H0)
    Lw0 <- Lw; Lw <- norm(Hs,'2') 
    gradW <- -X%*%t(H0) + W0%*%Hs
    W <- Wm - gradW/Lw
    W[W<0] <- 0
    # extract features: point multiplication, sum by row, sort
    ids <- order(apply(W*W,1,sum), decreasing = T)[1:k] 
    W[-ids, ] <- 0 # Row sparse
    
    # --------------------------------------------------------------------------
    # update H
    Ws <- t(W)%*%W + rho*matrix(1, ncol(W), ncol(W));
    Lh0 <- Lh; Lh <- norm(Ws + (v-rho)*diag(nrow(Ws)),'2') 
    gradH <- -t(W)%*%X + Ws%*%H0 + (v - rho)*H0
    H <- Hm - gradH/Lh;
    H[H<0] <- 0
    
    # --------------------------------------------------------------------------
    # calculate objective function value
    obj <- norm(X-W%*%H,'F')^2/2 + (v-rho)*norm(H,'F')^2/2 + rho*sum(apply(H,2,sum)^2)/2
    objs[iter] <- obj
    
    # --------------------------------------------------------------------------
    # stopping condition
    delta <- norm(W-W0,'F')/norm(W0,'F') + norm(H-H0,'F')/norm(H0,'F');
    if(iter>10 && delta <= palm_tol) break 
    
    # --------------------------------------------------------------------------
    # --- correction and extrapolation ---
    t <- (1+sqrt(1+4*t0^2))/2
    
    # --------------------------------------------------------------------------
    # restore to previous W,H, and cached quantities for nonincreasing objective
    if(obj>=obj0){
      Wm <- W0; Hm <- H0;
    }else{
      # extrapolation
      ww <- (t0-1)/t  # extrapolation weight
      wx <- min(c(ww,rw*sqrt(Lw0/Lw)))  # choose smaller one for convergence
      wy <- min(c(ww,rw*sqrt(Lw0/Lw)))
      
      Wm <- W + wx*(W-W0); Hm <- H + wy*(H-H0) # extrapolation
      
      W0 <- W; H0 <- H; t0 <- t; obj0 <- obj;
    }
  }
  return(list(W=W,H=H,rho=rho,objs=objs))
}