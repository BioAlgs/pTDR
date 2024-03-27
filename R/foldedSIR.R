######################################################################
## folded SIR
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


#' Folded Sliced Inverse Regression (SIR)
#'
#' Implements a folded version of Sliced Inverse Regression (SIR) for dimensionality reduction,
#' designed for bidimensional tensorial data. This function requires the data matrix `x` to
#' have dimensions that can be folded into a tensor of shape specified by `p`.
#'
#' @param y Numeric vector of response variables, with one value per observation.
#' @param x Numeric matrix of predictor variables, with each row representing an observation.
#' @param p Numeric vector of length 2 specifying the dimensions into which `x` can be folded.
#' @param nslice Integer specifying the number of slices to divide `y` into for SIR analysis.
#' @param d Numeric vector of length 2 specifying the dimensions for dimensionality reduction.
#' @param maxiter Maximum number of iterations for the algorithm (default = 200).
#' @param tol Tolerance level for algorithm convergence (default = 1e-8).
#'
#' @return A list containing:
#' \itemize{
#'   \item Various components related to the fitted model, including the number of iterations (`iter`) and convergence status (`conv`).
#' }
#'
#' @examples
#' set.seed(123)
#' y <- rnorm(100)
#' x <- matrix(rnorm(400), ncol = 4)
#' p <- c(2, 2)
#' nslice <- 10
#' d <- c(1, 1)
#'
#' result <- foldedSIR(y, x, p, nslice, d)
#' print(result)
#'
#' @export
foldedSIR = function(y, x, p, nslice, d, maxiter = 200, tol = 1e-8)
{
  if(length(p) != 2) stop("length(p) != 2.\n");
  if(length(d) != length(p)) stop("length(p) != length(d).\n");

  if(length(y) != NROW(x)) stop("length(y) != NROW(x).\n");
  if(prod(p) != NCOL(x)) stop("prod(p) != NCOL(x).\n");

  ## slice and calculate Sigma, Mmat

  info = SIR.slice(y, x, nslice);
  myprob = c(info, list(pL = p[1], pR = p[2], dL = d[1], dR = d[2]));

  ## initial values of a, b
  M0 = matrix(AB2eigen1st(myprob$Mmat, myprob$Smat)$vectors, nrow = myprob$pL);
  temp = svd(M0, nu = myprob$dL, nv = myprob$dR)
  a0 = matrix(temp$u, ncol = myprob$dL);
  b0 = matrix(temp$v, ncol = myprob$dR);
  boa0 = kronecker(b0, a0);
  V2V2 = t(boa0) %*% myprob$Smat %*% boa0;
  f0 = matrix(solve(V2V2, t(boa0) %*% t(myprob$slice.mean)), 
              nrow = prod(d), ncol = myprob$slice.info$nslice);

  ## start iteration
  newest = list(a = a0, b = b0, boa = boa0, f = f0, Qfun = 0.0);
  iter = 0; conv = tol + 1.0;
  while((iter < maxiter) && (conv > tol))
  {
    iter = iter + 1;
    oldest = newest;
    newest = foldedSIR.onestep(myprob, oldest);
    conv = abs(oldest$Qfun - newest$Qfun);
  }

  c(myprob, newest, list(iter = iter, conv = conv));
}

#' One Step of Folded SIR Iteration
#'
#' A helper function for `foldedSIR` that performs a single iteration step of the
#' folded SIR algorithm, updating the estimates of parameters.
#'
#' @param prob A list containing problem-specific parameters and data derived from the initial `foldedSIR` setup.
#' @param est A list containing the current estimates of the parameters.
#'
#' @return A list with updated parameter estimates.
#'
#' @noRd
foldedSIR.onestep = function(prob, est)
{
  newest = est;
 
  ## update a
  V2V2 = matrix(0.0, nrow = prob$pL*prob$dL, ncol = prob$pL*prob$dL);
  V2V1 = matrix(0.0, nrow = prob$pL*prob$dL, ncol = 1);
  for(i in 1:prob$slice.info$nslice)
  {
    bf = newest$b %*% matrix(newest$f[, i], nrow = prob$dR, byrow = TRUE);
    bfoI = kronecker(bf, diag(rep(1.0, prob$pL))); 
    V2V2 = V2V2 + prob$slice.prop[i] * t(bfoI) %*% prob$Smat %*% bfoI;
    V2V1 = V2V1 + prob$slice.prop[i] * t(bfoI) %*% prob$slice.mean[i, ];
  }
  a0 = matrix(solve(V2V2, V2V1), ncol = prob$dL)
  newest$a = matrix(colspace(a0), ncol = prob$dL);

  ## update b
  V2V2 = matrix(0.0, nrow = prob$pR*prob$dR, ncol = prob$pR*prob$dR);
  V2V1 = matrix(0.0, nrow = prob$pR*prob$dR, ncol = 1);
  for(i in 1:prob$slice.info$nslice)
  {
    af = newest$a %*% matrix(newest$f[, i], nrow = prob$dL);
    Ioaf = kronecker(diag(rep(1.0, prob$pR)), af);
    V2V2 = V2V2 + prob$slice.prop[i] * t(Ioaf) %*% prob$Smat %*% Ioaf;
    V2V1 = V2V1 + prob$slice.prop[i] * t(Ioaf) %*% prob$slice.mean[i, ];
  }
  b0 = matrix(solve(V2V2, V2V1), ncol = prob$dR);
  newest$b = matrix(colspace(b0), ncol = prob$dR);

  boa = kronecker(newest$b, newest$a);
  newest$boa = boa;

  ## update f
  V2V2 = t(boa) %*% prob$Smat %*% boa;
  newest$f = matrix(solve(V2V2, t(boa) %*% t(prob$slice.mean)),
                    nrow = prob$dL * prob$dR, 
                    ncol = prob$slice.info$nslice);

  ## update Qfun
  newest$Qfun = 0.0;
  for(i in 1:prob$slice.info$nslice)
  {
    tvec = prob$slice.mean[i, ] - boa %*% newest$f[, i];    
    newest$Qfun = newest$Qfun + prob$slice.prop[i] * t(tvec) %*% prob$Smat %*% tvec;
  }

  return(newest);
}

######################################################################
## THE END
######################################################################
