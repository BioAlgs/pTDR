######################################################################
## functions shared by tensorSIR and foldedSIR
######################################################################

require('dr');
dyn.load(paste("./inst/libs/d_eig", .Platform$dynlib.ext, sep = ""));
source('./R/locpoly.R')

######################################################################
## calculate cov(x), slice y, calculate var[E(x|y)] 
######################################################################

SIR.slice = function(y, x, nslice)
{
  xmean = colMeans(x);
  xcov = cov(x);

  sliceinfo = dr.slices(y, nslice);

  slicemean = matrix(0.0, nrow = sliceinfo$nslice, ncol = NCOL(x));
  sliceprop = rep(0.0, sliceinfo$nslice);

  Mmat = matrix(0, nrow = NCOL(x), ncol = NCOL(x));
  for(i in 1:sliceinfo$nslice)
  {
    index = (sliceinfo$slice.indicator == i);
    slicemean[i, ] = colMeans(x[index, ]) - xmean;
    sliceprop[i] = mean(index);
    Mmat = Mmat + sliceprop[i] * slicemean[i, ] %*% t(slicemean[i, ]);
  } 

  list(slice.info = sliceinfo, 
       slice.mean = slicemean, 
       slice.prop = sliceprop, 
       xmean = xmean, Smat = xcov, Mmat = Mmat);
}

######################################################################
## calculate eigenvalue and eigenvector of a symmetric matrix
## only lower triangle is used in calculation
######################################################################

Ceig = function(x, number = 1, tol = .Machine$double.eps)
{
  n = NROW(x); 
  if(NCOL(x) != n) stop("x should be a square matrix.\n");

  ans = .C("d_eig", as.integer(n), as.double(x), as.integer(number),
     as.double(tol), d = integer(1), 
     values = double(number), vectors = double(n*number));
  list(values  = ans$values[ans$d : 1],  
       vectors = matrix(ans$vectors, nrow = n)[, ans$d : 1])
}

######################################################################
## solve max_x (x^TAx) / (x^TBx)
## when B does not have full rank, we assume that x in span(B)
## assume that A and B are symmetric (low-triangle is used)
######################################################################

AB2eigen = function(Amat, Bmat, tol = 1e-7)
{
  Beig = eigen(Bmat, symmetric = TRUE);
  Brank = sum(abs(Beig$values) > tol);

  if(Brank > 1)
  {  
    Dmat = sweep(Beig$vectors[, 1:Brank], 2, sqrt(Beig$values[1:Brank]), "/");
  }
  else
  { ## Brank == 1
    Dmat = matrix(Beig$vectors[, 1] / sqrt(Beig$values[1]), ncol = 1);
  }

  Aeig = eigen(t(Dmat) %*% Amat %*% Dmat, symmetric = TRUE);
  ans = Dmat %*% Aeig$vectors;

  list(values = Aeig$values, vectors = ans, rank = Brank);
}

######################################################################
## similar to AB2eigen(), but only calculate 1st eigenvectors
######################################################################

AB2eigen1st = function(Amat, Bmat, tol = 1e-7)
{
  Beig = eigen(Bmat, symmetric = TRUE);
  Brank = sum(abs(Beig$values) > tol);

  if(Brank > 1)
  {  
    Dmat = sweep(Beig$vectors[, 1:Brank], 2, sqrt(Beig$values[1:Brank]), "/");
  }
  else
  { ## Brank == 1
    Dmat = matrix(Beig$vectors[, 1] / sqrt(Beig$values[1]), ncol = 1);
  }

  Aeig = Ceig(t(Dmat) %*% Amat %*% Dmat, number = 1);
  ans = Dmat %*% Aeig$vectors;

  list(values = Aeig$values, vectors = ans, rank = Brank);
}

######################################################################
## scale a vector to unit length
######################################################################

scale2unit = function(x)
{
  x / sqrt(sum(x * x));
}

######################################################################
## sqrt root of a square matrix, A^{1/2}
######################################################################

sqrtmat = function(A)
{
  Aeig = eigen(A, symmetric = TRUE);
  Aeig$vectors %*% diag(sqrt(Aeig$values)) %*% t(Aeig$vectors);
}

######################################################################
## inverse of the sqrt root of a square matrix, A^{-1/2}
######################################################################

sqrtmatinv = function(A)
{
  Aeig = eigen(A, symmetric = TRUE);
  Aeig$vectors %*% diag(1/sqrt(Aeig$values)) %*% t(Aeig$vectors);
}

######################################################################
## find the orthonormal basis of the column space of a matrix
######################################################################

colspace = function(x, tol = 1e-7)
{
  xqr = qr(x, tol = tol);
  qr.Q(xqr)[, 1:(xqr$rank)]
}

######################################################################
## find the orthonormal basis of the null space of a matrix
######################################################################

nullspace = function(x, tol = 1e-7)
{
  xqr = qr(x, tol = tol);

  if(xqr$rank == NROW(x))  ans = NULL
  else ans = qr.Q(xqr, complete = TRUE)[, (xqr$rank+1):NROW(x)];

  ans
}

######################################################################
## find the projection matrix of a matrix
######################################################################

projmat = function(x, tol = 1e-7)
{
  xqr = qr(x, tol = tol);
  Q = qr.Q(xqr)[, 1:(xqr$rank)];
  Q %*% t(Q);
}

######################################################################
## distance between the column spaces of two matrices
## Suppose A and B are orthonormal matrix,
## the distance is the average of eigenvalue (r)
##
## return NA if matA or matB contains NA
##           if matA or matB is almost zero
######################################################################

dist.colspace = function(matA, matB)
{
  if(any(is.na(matA)) || any(is.na(matB)))
  {
    ans = NA
  } 
  else if((max(abs(matA)) < 1e-15) || (max(abs(matB)) < 1e-15))
  {
    ans = NA
  }
  else
  {
    A.orth = qr.Q(qr(matA))
    B.orth = qr.Q(qr(matB))
    BAAB = t(B.orth) %*% A.orth %*% t(A.orth) %*% B.orth
    BAAB.eig = eigen(BAAB, only.values = T)$values
    ans = 1 - sqrt(mean(BAAB.eig))
  }

  ans
}

######################################################################
## calculate Qk = I - S %*% M %*% solve(t(M) %*% S %*% M) %*% t(M)
######################################################################

cal.Qkmat = function(Sigma, Mkmat)
{
  SM = Sigma %*% Mkmat;
  ans = - SM %*% solve(t(Mkmat) %*% SM, t(Mkmat));
  diag(ans) = diag(ans) + 1.0;
  ans;
}

######################################################################
## cross-validation for given B
######################################################################

CV4B = function(y, x, Bmat, h, m = 10, group = NULL)
{
  n = length(y)
  group = if(is.null(group)) sample((1:n)%%m) + 1 else as.numeric(factor(group));

  xB = x %*% Bmat; CVscore = numeric(length(unique(group)));
  for(i in 1:m)
  {
    index = (group == i);
    x0 = xB[index, ]; y0 = y[index];
    x1 = xB[!index, ]; y1 = y[!index];
    yhat = Clocpoly(y1, x1, h, x0);
    CVscore[i] = mean((y0 - yhat)^2);
  }

  list(group = group, CVscore = CVscore, CV = mean(CVscore), CVstd = sd(CVscore));
}

######################################################################
## THE END
######################################################################
