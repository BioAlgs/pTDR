#' Tensor Reduction and Estimation Package
#' @title Tensor Reduction and Estimation for TensorSIR
#' @description
#' This package provides tools for performing tensor reduction and estimation, 
#' leveraging both R and C  for computational efficiency. Tensor operations can 
#' be highly computational intensive, especially with large datasets commonly 
#' found in machine learning, data mining, and high-dimensional data analysis.
#' 
#' @name TensorSIR
#'
#' @examples
#' ######################################################################
#' ## example
#' ######################################################################
#' rm(list = ls());
#' a1 = c(1, 1, 0, 0, 0);
#' a2 = c(0, 0, 0, 1, 1);
#' amat = cbind(a1, a2);
#' b1 = c(1, 1, 1, rep(0, 6));
#' b2 = c(rep(0, 6), 1, 1, 1);
#' bmat = cbind(b1, b2);
#' 
#' x = matrix(rnorm(400 * 45), ncol = 45);
#' boa = kronecker(bmat, amat)
#' y = (x %*% boa[, 1]) / (2 + (3 + x %*% boa[, 4])^2) + rnorm(400, sd = 0.5);
#' 
#' p = c(5, 9);
#' d = c(1, 1);
#' nslice = 10;
#' 
#' ans1 = foldedSIR(y, x, c(5, 9), 10, c(2, 2));
#' ans2 = pTDR(y, x, c(5, 9), 10, 4);
#' 
#' dist.colspace(boa, ans1$boa)
#' dist.colspace(boa, ans2$direction)
#' 
#' #####################################################################
#' ## compare the new old with the old one
#' #####################################################################
#' 
#' rm(list = ls());
#' # source('./R/foldedSIR.R');
#' # source('./R/shared-funs.R');
#' # source('./R/foldedSIR.R');
#' 
#' a1 = c(1, 1, 0, 0, 0);
#' a2 = c(0, 0, 0, 1, 1);
#' amat = cbind(a1, a2);
#' b1 = c(1, 1, 1, rep(0, 6));
#' b2 = c(rep(0, 6), 1, 1, 1);
#' bmat = cbind(b1, b2);
#' 
#' x = matrix(rnorm(400 * 45), ncol = 45);
#' boa = kronecker(bmat, amat)
#' y = (x %*% boa[, 1]) / (2 + (3 + x %*% boa[, 4])^2) + rnorm(400, sd = 0.5);
#' 
#' ans0 = foldedSIR(y, x, c(5, 9), 10, c(1, 1));
#' ans1 = foldedSIR(y, x, c(5, 9), 10, c(1, 1));
#' 
#' dist.colspace(ans0$a, ans1$a)
#' dist.colspace(ans0$b, ans1$b)
#' 
#' boa0 = kronecker(ans0$b, ans0$a);
#' boa1 = kronecker(ans1$b, ans1$a);
#' dist.colspace(boa0, boa1);
#' 
#' ######################################################################
#' ## testing
#' ######################################################################
#' 
#' rm(list = ls());
#' 
#' a1 = c(1, 1, 0, 0, 0);
#' a2 = c(0, 0, 0, 1, 1);
#' amat = cbind(a1, a2);
#' b1 = c(1, 1, 1, rep(0, 6));
#' b2 = c(rep(0, 6), 1, 1, 1);
#' bmat = cbind(b1, b2);
#' 
#' x = matrix(rnorm(400 * 45), ncol = 45);
#' boa = kronecker(bmat, amat)
#' y = (x %*% boa[, 1]) / (2 + (3 + x %*% boa[, 4])^2) + rnorm(400, sd = 0.5);
#' 
#' ans = testH0(y, x, 10);
#' ans$p.value
#' wchisq.tail(ans$stat, ans$z);
#' 
#' ######################################################################
#' ## cross-validation
#' ######################################################################
#' 
#' rm(list = ls());
#' 
#' a1 = c(1, 1, 0, 0, 0);
#' a2 = c(0, 0, 0, 1, 1);
#' amat = cbind(a1, a2);
#' b1 = c(1, 1, 1, rep(0, 6));
#' b2 = c(rep(0, 6), 1, 1, 1);
#' bmat = cbind(b1, b2);
#' 
#' x = matrix(rnorm(400 * 45), ncol = 45);
#' boa = kronecker(bmat, amat)
#' y = (x %*% boa[, 1]) / (2 + (3 + x %*% boa[, 4])^2) + rnorm(400, sd = 0.5);
#' 
#' ans = CV4B(y, x, boa[, c(1, 4)], 0.5)
#' 
NULL