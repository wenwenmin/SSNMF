NMF = function(M, X0, Y0, tol=3e-3, maxit=500){
  
  # NMF: nonnegative matrix factorization by block-coordinate update
  # minimize 0.5*||M - X*Y||_F^2
  # subject to X>=0, Y>=0
  
  # input:
  #       M: input nonnegative matrix (p times n)
  #       r: rank (X has r columns and Y has r rows)

  # Reference: 
  # Y. Xu and W. Yin, A block coordinate descent method for regularized 
  # multiconvex optimization with applications to nonnegative tensor factorization 
  # and completion, SIAM Journal on Imaging Sciences, 6(3), 1758-1789, 2013
  # More information can be found at:
  # https://www.math.ucla.edu/~wotaoyin/papers/bcu/nmf/nmf.html
  
  rw = 1;
  ## Data preprocessing and initialization
  m = nrow(M)
  n = ncol(M)
  
  Mnrm = norm(M,'f')
  
  X0 = X0/norm(X0,'f')*sqrt(Mnrm)
  Y0 = Y0/norm(Y0,'f')*sqrt(Mnrm)
  
  Xm = X0 
  Ym = Y0 
  Yt = t(Y0)
  Ys = Y0%*%Yt
  MYt = M%*%Yt;
  
  obj0 = 0.5*Mnrm*Mnrm;
  
  nstall = 0
  t0 = 1
  Lx = 1
  Ly = 1
  
  # cat("nonnegative matrix factorization by block-coordinate update\n")
  objs <- c() # Objective function values
  for(k in 1:maxit){
    # Gx, Gy: gradients with respect to X, Y
    # X, Y: new updates
    # Xm, Ym: extrapolations of X, Y
    # Lx, Lx0: current and previous Lipschitz bounds used in X-update
    # Ly, Ly0: current and previous Lipschitz bounds used in Y-update
    # obj, obj0: current and previous objective values
    # Xs = X'*X, Ys = Y*Y'
    
    # --- X-update ---
    Lx0 = Lx; Lx  = norm(Ys,"f");  # save and update Lipschitz bound for X
    Gx = Xm%*%Ys - MYt; # gradient at X=Xm
    # X = max(0, Xm - Gx/Lx)
    X = Xm - Gx/Lx; X[X<0] = 0
    Xt = t(X); Xs = Xt%*%X; # cache useful computation
    
    # --- Y-update ---
    Ly0 = Ly; Ly = norm(Xs, "f");  # save and update Lipschitz bound for Y
    Gy = Xs%*%Ym - Xt%*%M;    # gradient at Y=Ym
    # Y = max(0, Ym - Gy/Ly);
    Y = Ym - Gy/Ly; Y[Y<0] = 0
    Yt = t(Y); Ys = Y%*%Yt; MYt = M%*%Yt;  # cache useful computation
    
    # --- diagnostics, reporting, stopping checks ---
    # obj = 0.5*(sum(sum(Xs.*Ys)) - 2*sum(sum(X.*MYt)) + Mnrm*Mnrm);
    obj = 0.5*norm(M-X%*%Y,"f")^2  
    
    objs[k] = obj
    # stall and stopping checks
    
    # stopping condition
    delta <- norm(X-X0,'F')/norm(X0,'F') + norm(Y-Y0,'F')/norm(Y0,'F');
    if(k > 20 && delta <= tol) break
    
    # --- correction and extrapolation ---
    t = (1+sqrt(1+4*t0^2))/2
    
    if(obj>=obj0){
      # restore to previous X,Y, and cached quantities for nonincreasing objective
      Xm = X0; Ym = Y0;
      Yt = t(Y0); Ys = Y0%*%Yt; MYt = M%*%Yt;
    }else{
      # extrapolation
      w = (t0-1)/t  # extrapolation weight
      wx = min(c(w,rw*sqrt(Lx0/Lx)))  # choose smaller one for convergence
      wy = min(c(w,rw*sqrt(Ly0/Ly)))
      Xm = X + wx*(X-X0); Ym = Y + wy*(Y-Y0) # extrapolation
      X0 = X; Y0 = Y; t0 = t; obj0 = obj;
    }
  }
  return(list(W=X, H=Y, objs=objs))
}