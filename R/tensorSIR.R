######################################################################
## tensor SIR
## 02-05-2015
## Peng Zeng
######################################################################

######################################################################
## The response is a scalar and the predictor is an m-by-p matrix.
## Assume that the sample size is n
## Input: 
##    y --- n-by-1 vector
##    x --- n-by-(pL-times-pR) matrix
##    p --- c(pL, pR)
##    nslice --- number of slices
##    d --- c(dL, dR)
######################################################################

tensorSIR = function(y, x, p, nslice, ndir = 1, 
     maxiter = 200, tol = 1e-8, eps = 1e-8)
{
  if((length(p) < 1) || (length(p) > 3)) stop("length(p) != 1, 2, 3.\n");
  if(length(y) != NROW(x)) stop("length(y) != NROW(x).\n");
  if(prod(p) != NCOL(x)) stop("prod(p) != NCOL(x).\n");

  info = SIR.slice(y, x, nslice);
  myprob = c(info, list(p = p, ndir = ndir, 
                        maxiter = maxiter, tol = tol, eps = eps));

  ans = switch(length(p),
               tensor1SIR(myprob),
               tensor2SIR(myprob),
               tensor3SIR(myprob));

  c(list(p = p, ndir = ndir), info, ans);
}

######################################################################
## tensorSIR for vector-valued predictors
######################################################################

tensor1SIR = function(prob)
{
  if(length(prob$p) != 1) stop("length(p) != 1.\n");

  ans = AB2eigen(prob$Mmat, prob$Smat);
  ans$vectors = apply(ans$vectors, 2, scale2unit);

  list(direction = ans$vectors[, 1:prob$ndir],
       Qvalue = ans$values[1:prob$ndir]);
}

######################################################################
## tensorSIR for matrix-valued predictors: naive algorithm
######################################################################

tensor2SIR = function(prob)
{
  if(length(prob$p) != 2) stop("length(p) != 2.\n");

  amat = NULL; bmat = NULL; Qvalue = NULL; boa = NULL;
  Qk = diag(rep(1.0, prod(prob$p)));

  subprob = c(prob, list(pL = prob$p[1], pR = prob$p[2], 
                         a0 = NULL, b0 = NULL));

  for(i in 1:prob$ndir)
  {
    subprob$Mmat = Qk %*% prob$Mmat %*% t(Qk);   
    subprob$Smat = Qk %*% prob$Smat %*% t(Qk);  
    
    istep = tensor2SIR.pair.1st(subprob);

    amat = cbind(amat, istep$a);
    bmat = cbind(bmat, istep$b); 
    boa  = cbind(boa,  istep$boa);
    Qvalue = c(Qvalue, istep$Qfun);

    ## find the projection matrix for the next pair of directions 
    Qk = cal.Qkmat(prob$Smat, boa);
  }

  list(direction = boa, Qvalue = Qvalue, a = amat, b = bmat);
}

######################################################################
## extract one pair of directions 
######################################################################

tensor2SIR.pair.1st = function(prob)
{
  newest = list(pL = prob$pL, pR = prob$pR, 
                a = prob$a0, b = prob$b0, boa = NULL, Qfun = 0);

  if(max(abs(prob$Mmat)) < prob$eps) return( newest );
  ## stop if Mmat == 0.

  if(is.null(prob$a0) || is.null(prob$b0))
  {
    M0 = matrix(AB2eigen1st(prob$Mmat, prob$Smat)$vectors, nrow = prob$pL);
    temp = svd(M0, nu = 1, nv = 1)
    a0 = matrix(temp$u, ncol = 1);
    b0 = matrix(temp$v, ncol = 1);
  }
  newest$a = scale2unit( a0 );
  newest$b = scale2unit( b0 );

  iter = 0; conv = prob$tol + 1.0;
  while((iter < prob$maxiter) && (conv > prob$tol))
  {
    iter = iter + 1;
    oldest = newest;
    newest = tensor2SIR.onestep(prob$Mmat, prob$Smat, oldest);
    conv = abs(oldest$Qfun - newest$Qfun);
  }

  if((iter >= prob$maxiter) && (conv > prob$tol))
    warning("Reach maxiter, increase maxiter and/or increase tol.\n");

  c(newest, list(iter = iter, conv = conv));
}

######################################################################
## one step. update a, b, and Q
######################################################################

tensor2SIR.onestep = function(Mmat, Smat, est, tol = 1e-7)
{
  ans = est;

  ## update b when pR > 1

  if(ans$pR > 1)
  {
    Ia = kronecker(diag(rep(1, ans$pR)), ans$a);
    Ma = t(Ia) %*% Mmat %*% Ia;
    Sa = t(Ia) %*% Smat %*% Ia;
    atmp = AB2eigen1st(Ma, Sa);

    ans$b = scale2unit(atmp$vector[, 1]);
    ans$Qfun = atmp$value[1];
  } 
  else ans$b = matrix(1.0, nrow = 1, ncol = 1);

  ## update a when pL > 1

  if(ans$pL > 1)
  {
    bI = kronecker(ans$b, diag(rep(1, ans$pL)));
    Mb = t(bI) %*% Mmat %*% bI;
    Sb = t(bI) %*% Smat %*% bI;
    btmp = AB2eigen1st(Mb, Sb);

    ans$a = scale2unit(btmp$vector[, 1]);
    ans$Qfun = btmp$value[1];
  } 
  else ans$a = matrix(1.0, nrow = 1, ncol = 1);
  
  ## update boa and Qfun

  ans$boa = kronecker(ans$b, ans$a);

  baSba = drop( t(ans$boa) %*% Smat %*% ans$boa );
  if(abs(baSba) < tol) 
  {
    ans$Qfun = 0.0 
  }
  else 
  {
    ans$Qfun = drop( (t(ans$boa) %*% Mmat %*% ans$boa) / baSba );
  }

  ## output

  ans;
}

######################################################################
## THE END
######################################################################
