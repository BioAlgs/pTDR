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

#' Tensor-based Sliced Inverse Regression (SIR)
#'
#' Implements a tensor-based approach to Sliced Inverse Regression (SIR) for 
#' dimensionality reduction in supervised learning scenarios. This function 
#' can handle tensor inputs of up to 3 dimensions, specified by the `p` parameter.
#'
#' @param y Numeric vector of response variables, indicating the outcome for each observation.
#' @param x Numeric matrix of predictor variables, where each row is an observation 
#'   and columns correspond to variables, possibly representing a flattened tensor.
#' @param p Numeric vector indicating the dimensions of the tensor representation 
#'   for the predictor variables. The product of `p` should equal the number of columns in `x`.
#' @param nslice Integer specifying the number of slices for the SIR analysis.
#' @param ndir Integer indicating the number of directions to estimate (default is 1).
#' @param maxiter Maximum number of iterations for the algorithm (default is 200).
#' @param tol Tolerance level for convergence in the optimization algorithm (default is 1e-8).
#' @param eps A small numerical value to stabilize computations (default is 1e-8).
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item p: The dimensions of the tensor as specified by the input.
#'   \item ndir: The number of directions estimated.
#'   \item info: Information from the initial SIR analysis, including slice information.
#'   \item ans: Results from the tensor-based SIR analysis, varying based on the dimensionality of `p`.
#' }
#'
#' @examples
#' set.seed(123)
#' y <- rnorm(100)
#' x <- matrix(rnorm(300), ncol = 3)
#' p <- c(1, 3)
#' nslice <- 10
#' 
#' result <- pTDR(y, x, p, nslice)
#' print(result)
#'
#' @export
pTDR = function(y, x, p, nslice, ndir = 1, 
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
## pTDR for vector-valued predictors
######################################################################

#' Tensor-based SIR for Vector-valued Predictors
#'
#' Performs tensor-based Sliced Inverse Regression (SIR) for cases where predictors
#' are vector-valued, utilizing a single dimension (`p` = 1).
#'
#' @param prob A list containing problem-specific parameters, including the
#'   moment matrices `Mmat` and `Smat`, dimension `p`, and the number of directions `ndir`.
#'
#' @return A list containing the estimated directions (`direction`) and their 
#'   corresponding eigenvalues (`Qvalue`).
#'
#' @examples
#' # For usage, see the example for `foldedSIR`.
#'
#' @noRd
tensor1SIR = function(prob)
{
  if(length(prob$p) != 1) stop("length(p) != 1.\n");

  ans = AB2eigen(prob$Mmat, prob$Smat);
  ans$vectors = apply(ans$vectors, 2, scale2unit);

  list(direction = ans$vectors[, 1:prob$ndir],
       Qvalue = ans$values[1:prob$ndir]);
}

######################################################################
## pTDR for matrix-valued predictors: naive algorithm
######################################################################

#' Tensor-based SIR for Matrix-valued Predictors Using Naive Algorithm
#'
#' Extends tensor-based Sliced Inverse Regression (SIR) to handle matrix-valued predictors
#' by performing the analysis in a tensor framework, specifically for two dimensions (`p` = 2).
#'
#' @param prob A list with problem-specific parameters and moment matrices.
#'
#' @return A list containing the tensor-based directions (`direction`), their Q-values (`Qvalue`),
#'   and matrices `a` and `b` representing the tensor decomposition.
#'
#' @examples
#' # For usage, see the example for `foldedSIR`.
#'
#' @noRd
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

#' Extract One Pair of Directions for Tensor2SIR
#'
#' A helper function for `tensor2SIR` to extract a single pair of directions in the 
#' tensor-based SIR analysis for matrix-valued predictors.
#'
#' @param prob A list containing the modified problem-specific parameters for the current iteration.
#'
#' @return A list with the newly estimated parameters, including the number of iterations (`iter`)
#'   and convergence status (`conv`).
#'
#' @noRd
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

#' One Step Update for Tensor2SIR Pair Extraction
#'
#' Performs a single step update of parameters `a` and `b` in the tensor2SIR algorithm,
#' aiming to refine the direction estimates for matrix-valued predictors.
#'
#' @param Mmat Moment matrix M associated with the current slice.
#' @param Smat Covariance matrix S of the predictors.
#' @param est A list containing the current estimates of `a` and `b`.
#' @param tol Tolerance for numerical stability checks (default = 1e-7).
#'
#' @return A list with updated estimates of `a`, `b`, their outer product (`boa`), and `Qfun`.
#'
#' @noRd
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
