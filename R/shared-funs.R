######################################################################
## functions shared by tensorSIR and foldedSIR
######################################################################

require('dr');
dyn.load(paste("./inst/libs/d_eig", .Platform$dynlib.ext, sep = ""));
# source('./R/locpoly.R')

######################################################################
## calculate cov(x), slice y, calculate var[E(x|y)] 
######################################################################

#' Slice Data for SIR Analysis
#'
#' Prepares sliced data for Sliced Inverse Regression (SIR) analysis, including computing 
#' means and proportions for each slice, as well as the necessary matrices (`Mmat` and `Smat`).
#'
#' @param y Numeric vector of response variables.
#' @param x Numeric matrix of predictor variables.
#' @param nslice Integer specifying the number of slices.
#'
#' @return A list containing slicing information, mean and proportion of each slice, 
#'   the mean of `x`, covariance matrix `Smat`, and matrix `Mmat` for SIR analysis.
#'
#' @noRd
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

#' Calculate Eigenvalues and Eigenvectors of Symmetric Matrix
#'
#' Computes the eigenvalues and eigenvectors of a symmetric matrix using its lower triangle.
#' Optimized by a C interface for efficiency.
#'
#' @param x Symmetric square matrix.
#' @param number Integer, the number of eigenvalues and vectors to compute.
#' @param tol Numeric, tolerance for determining positive values.
#'
#' @return A list with eigenvalues (`values`) and eigenvectors (`vectors`).
#'
#' @noRd
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

#' Solve Maximized Ratio of Quadratic Forms
#'
#' Solves the maximization problem (x^TAx) / (x^TBx) given symmetric matrices A and B,
#' assuming x in the span of B and utilizing eigen decomposition.
#'
#' @param Amat Numeric matrix, symmetric.
#' @param Bmat Numeric matrix, symmetric, not necessarily full rank.
#' @param tol Numeric, tolerance for rank determination of B.
#'
#' @return A list containing the eigenvalues (`values`), the eigenvectors (`vectors`),
#'   and the rank of B (`rank`).
#'
#' @noRd
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

#' Calculate First Eigenvector for Maximized Ratio of Quadratic Forms
#'
#' Specifically targets the first eigenvector in the maximization problem (x^TAx) / (x^TBx),
#' for symmetric matrices A and B, considering x in the span of B.
#'
#' @param Amat Numeric matrix, symmetric.
#' @param Bmat Numeric matrix, symmetric, not necessarily full rank.
#' @param tol Numeric, tolerance for rank determination of B.
#'
#' @return A list with the first eigenvalue (`values`), the corresponding eigenvector (`vectors`),
#'   and the rank of B (`rank`).
#'
#' @noRd
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

#' Scale Vector to Unit Length
#'
#' Scales a numeric vector so that its length is 1, based on its Euclidean norm.
#'
#' @param x Numeric vector to be scaled.
#'
#' @return Numeric vector scaled to unit length.
#'
#' @noRd
scale2unit = function(x)
{
  x / sqrt(sum(x * x));
}

######################################################################
## sqrt root of a square matrix, A^{1/2}
######################################################################

#' Square Root of a Symmetric Matrix
#'
#' Computes the square root of a symmetric matrix.
#'
#' @param A Symmetric square matrix.
#'
#' @return Square root of matrix A.
#'
#' @noRd
sqrtmat = function(A)
{
  Aeig = eigen(A, symmetric = TRUE);
  Aeig$vectors %*% diag(sqrt(Aeig$values)) %*% t(Aeig$vectors);
}

######################################################################
## inverse of the sqrt root of a square matrix, A^{-1/2}
######################################################################

#' Inverse Square Root of a Symmetric Matrix
#'
#' Computes the inverse square root of a symmetric matrix.
#'
#' @param A Symmetric square matrix.
#'
#' @return Inverse square root of matrix A.
#'
#' @noRd
sqrtmatinv = function(A)
{
  Aeig = eigen(A, symmetric = TRUE);
  Aeig$vectors %*% diag(1/sqrt(Aeig$values)) %*% t(Aeig$vectors);
}

######################################################################
## find the orthonormal basis of the column space of a matrix
######################################################################

#' Orthonormal Basis of the Column Space
#'
#' Finds the orthonormal basis of the column space of a matrix.
#'
#' @param x Matrix for which the column space is to be computed.
#' @param tol Tolerance for determining numerical rank.
#'
#' @return Orthonormal basis of the column space of x.
#'
#' @noRd
colspace = function(x, tol = 1e-7)
{
  xqr = qr(x, tol = tol);
  qr.Q(xqr)[, 1:(xqr$rank)]
}

######################################################################
## find the orthonormal basis of the null space of a matrix
######################################################################

#' Orthonormal Basis of the Null Space
#'
#' Finds the orthonormal basis of the null space of a matrix.
#'
#' @param x Matrix for which the null space is to be computed.
#' @param tol Tolerance for determining numerical rank.
#'
#' @return Orthonormal basis of the null space of x.
#'
#' @noRd
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

#' Projection Matrix of a Matrix
#'
#' Computes the projection matrix associated with a given matrix.
#'
#' @param x Matrix for which the projection matrix is to be computed.
#' @param tol Tolerance for determining numerical rank.
#'
#' @return Projection matrix of x.
#'
#' @noRd
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

#' Distance Between Column Spaces
#'
#' Calculates the distance between the column spaces of two orthonormal matrices.
#'
#' @param matA First orthonormal matrix.
#' @param matB Second orthonormal matrix.
#'
#' @return Distance between the column spaces of matA and matB.
#'
#' @noRd
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

#' Calculate Qk Matrix
#'
#' Computes the Qk matrix for dimension reduction, given a covariance matrix and a matrix M.
#'
#' @param Sigma Covariance matrix.
#' @param Mkmat Matrix M used in the calculation.
#'
#' @return Qk matrix.
#'
#' @noRd
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

#' Cross-Validation for Bandwidth Selection
#'
#' Performs m-fold cross-validation to select the optimal bandwidth for local polynomial regression. 
#' The function divides the data into `m` groups, either randomly or based on a provided grouping, 
#' and calculates the cross-validation score for the given bandwidth `h`.
#'
#' @param y Numeric vector of response variables.
#' @param x Numeric matrix of predictor variables.
#' @param Bmat Transformation matrix applied to predictors before regression.
#' @param h Bandwidth parameter for local polynomial regression.
#' @param m Number of groups for cross-validation (default is 10).
#' @param group Optional vector specifying the grouping of observations. 
#' If `NULL`, observations are randomly divided into groups.
#'
#' @return A list containing:
#' \itemize{
#'   \item{group}{Vector indicating the group assignment of each observation.}
#'   \item{CVscore}{Numeric vector of cross-validation scores for each group.}
#'   \item{CV}{Mean of the cross-validation scores.}
#'   \item{CVstd}{Standard deviation of the cross-validation scores.}
#' }
#'
#' @examples
#' # Simulated data
#' set.seed(123)
#' y <- rnorm(100)
#' x <- matrix(rnorm(200), ncol = 2)
#' Bmat <- diag(2)
#' h <- 1
#'
#' # Perform cross-validation
#' cv_results <- CV4B(y, x, Bmat, h)
#' print(cv_results)
#'
#' @export
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
