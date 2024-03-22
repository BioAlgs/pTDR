######################################################################
## SIR for matrix-valued predictors. one direction on left/right
## last updated: 01-09-2013 by Peng Zeng
######################################################################

######################################################################
## The response is a scalar and the predictor is an m-by-p matrix.
## Assume that the sample size is n
## Input: 
##    y --- n-by-1 vector
##    x --- n-by-(pL-times-pR) matrix
##    nslice --- number of slices
######################################################################

######################################################################
## the algorithm proposed in Li, Kim, and Altman (2010)
## only works when d = c(1, 1)
######################################################################

foldedSIR.old = function(y, x, p, nslice, d, maxiter = 200, tol = 1e-8)
{
  if(length(p) != 2) stop("length(p) != 2.\n");
  if(length(d) != length(p)) stop("length(p) != length(d).\n");

  if(length(y) != NROW(x)) stop("length(y) != NROW(x).\n");
  if(prod(p) != NCOL(x)) stop("prod(p) != NCOL(x).\n");

  if(any(d != c(1, 1))) stop("Currently only work when d = c(1, 1).\n");

  ## assemble inputs for foldedSIR
  info = SIR.slice(y, x, nslice);
  prob = c(list(type = "solve_foldedSIR",
              p = p, d = d, Smat = info$Smat, Mmat = info$Mmat,
              maxiter = maxiter, tol = tol), info);

  subprob = c(prob, list(pL = prob$p[1], pR = prob$p[2], 
                         dL = prob$d[1], dR = prob$d[2])); 
  subprob$type = "solve_foldedSIR211";

  ans = foldedSIR211(subprob);
  c(info, list(p = prob$p, d = prob$d), ans);
}


foldedSIR211 = function(prob)
{
  if(prob$type != "solve_foldedSIR211") stop("incorrect type.\n");
  if(length(prob$p) != 2) stop("length(p) != 2.\n");
  if(any(prob$d != c(1, 1))) stop("Currently only work when d = c(1, 1).\n");

  Sigma2    = sqrtmat(prob$Smat);
  Sigma2inv = sqrtmatinv(prob$Smat);

  ## find the initial value of a and b
  Mvector = Ceig(prob$Mmat, number = 1)$vectors;
  temp = svd(matrix(Mvector, nrow = prob$pL, ncol = prob$pR), nu = 1, nv = 1)
  a0 = scale2unit( matrix(temp$u, ncol = 1) );
  b0 = scale2unit( matrix(temp$v, ncol = 1) );

  f0 = numeric(prob$slice.info$nslices);
  ba = kronecker(b0, a0);
  V2V2 = drop( t(ba) %*% prob$Smat %*% ba );
  for(i in 1:prob$slice.info$nslices)  
  {
    V1V2 = t(a0) %*% matrix(prob$slice.mean[i, ], nrow = prob$pL) %*% b0;
    f0[i] = V1V2 / V2V2;
  }

  newest = c(prob, list(Sigma2 = Sigma2, Sigma2inv = Sigma2inv,
             a = a0, b = b0, boa = ba, f = f0, Qfun = 0));

  iter = 0; conv = prob$tol + 1;
  while((iter < prob$maxiter) && (conv > prob$tol))
  {
    iter = iter + 1;
    oldest = newest;
    newest = foldedSIR211.onestep(newest);
    conv = abs(oldest$Qfun - newest$Qfun);
  }

  c(newest, list(iter = iter, conv = conv, direction = newest$boa));
}


foldedSIR211.onestep = function(oldest)
{
  est = oldest;
  pfX = sweep(est$slice.mean, 1, est$slice.prop * est$f, "*");

  ## update b.

  Ia = kronecker(diag(rep(1, est$pR)), est$a);
  V2V2 = (t(Ia) %*% est$Smat %*% Ia) * sum(est$slice.prop * est$f^2);
  V1V2 = drop( t(est$a) %*% matrix(colSums(pfX), nrow = est$pL) );
  est$b = scale2unit( solve(V2V2, V1V2) );

  ## update a.

  bI = kronecker(est$b, diag(rep(1, est$pL)));
  V2V2 = (t(bI) %*% est$Smat %*% bI) * sum(est$slice.prop * est$f^2);
  V1V2 = drop( matrix(colSums(pfX), nrow = est$pL) %*% est$b );
  est$a = scale2unit( solve(V2V2, V1V2) );

  ## update f.

  est$boa = kronecker(est$b, est$a);
  V2V2 = drop( t(est$boa) %*% est$Smat %*% est$boa );
  for(i in 1:est$slice.info$nslices)  
  {
    V1V2 = t(est$a) %*% matrix(est$slice.mean[i, ], nrow = est$pL) %*% est$b;
    est$f[i] = V1V2 / V2V2;
  }

  ## update Qfun

  est$Qfun = 0;
  for(i in 1:est$slice.info$nslices)
  {
    tvec = est$Sigma2inv %*% est$slice.mean[i, ] - (est$Sigma2 %*% est$boa) * est$f[i];
    est$Qfun = est$Qfun + est$slice.prop[i] * sum(tvec^2);
  }

  ## output
  
  est;
}

######################################################################
## sqrt root of a square matrix, A^{1/2}
######################################################################

sqrtmat = function(A, tol = 1e-8)
{
  Aeig = eigen(A, symmetric = TRUE);
  Apos = (Aeig$values > tol)
  Aeig$vectors[, Apos] %*% diag(sqrt(Aeig$values[Apos])) %*% t(Aeig$vectors[, Apos]);
}

######################################################################
## inverse of the sqrt root of a square matrix, A^{-1/2}
######################################################################

sqrtmatinv = function(A, tol = 1e-8)
{
  Aeig = eigen(A, symmetric = TRUE);
  Apos = (Aeig$values > tol)
  Aeig$vectors[, Apos] %*% diag(1/sqrt(Aeig$values[Apos])) %*% t(Aeig$vectors[, Apos]);
}

######################################################################
## THE END
######################################################################
