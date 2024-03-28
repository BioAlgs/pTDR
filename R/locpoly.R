#' @useDynLib locpoly
#' @importFrom dr dr.slices
#' @export
dr::dr.slices

######################################################################
## local polynomial smoothing
######################################################################

dyn.load(paste("inst/libs/x64/locpoly", .Platform$dynlib.ext, sep = ""));


#' Local Polynomial Regression Using C Interface
#'
#' Performs local polynomial regression via a C function interface. It computes the
#' estimated values of a dependent variable `yin` based on the predictors `xin`, 
#' for new observations `xout`, using bandwidth `h`.
#'
#' @param yin Numeric vector of dependent variable values.
#' @param xin Numeric matrix of predictor variables.
#' @param h Numeric, the bandwidth for the local polynomial regression.
#' @param xout Numeric matrix specifying the points at which predictions are to be made.
#'
#' @return Numeric vector of predicted values corresponding to `xout`.
#'
#' @noRd
Clocpoly = function(yin, xin, h, xout)
{
  n = length(yin); p = NCOL(xin); m = NROW(xout);
  if(n != NROW(xin))  stop('NROW(xin) != length(yin).\n');
  if(p != NCOL(xout)) stop('NCOL(xin) != NCOL(xout).\n');

  .C('locpoly', as.integer(n), as.integer(p), as.integer(m), 
      as.double(yin), as.double(xin), as.double(h),
      as.double(xout), yhat = double(m),
     PACKAGE = "locpoly")$yhat;
}

#' Local Polynomial Regression in R
#'
#' An R implementation of local polynomial regression. It estimates values of the
#' dependent variable `yin` based on predictors `xin` for new observations `xout`,
#' using a specified bandwidth `h`.
#'
#' @param yin Numeric vector of dependent variable values.
#' @param xin Numeric matrix of predictor variables.
#' @param h Numeric, the bandwidth for the local polynomial regression.
#' @param xout Numeric matrix specifying the points at which predictions are to be made.
#'
#' @return Numeric vector of predicted values corresponding to `xout`.
#'
#' @noRd
Rlocpoly = function(yin, xin, h, xout)
{
  n = length(yin); p = NCOL(xin); m = NROW(xout);
  if(n != NROW(xin))  stop('NROW(xin) != length(yin).\n');
  if(p != NCOL(xout)) stop('NCOL(xin) != NCOL(xout).\n');

  yhat = numeric(m);
  for(k in 1:m)
  {
    xij = scale(xin, center = xout[k, ], scale = FALSE);
    wvec = exp((-0.5 / h / h) * rowSums(xij^2));
    wvec = wvec / sum(wvec);

    myfit = lm.wfit(cbind(1, xij), yin, wvec);
    yhat[k] = myfit$coefficients[1];
  }

  yhat;
}

######################################################################
## rule-of-thumb bandwidth
######################################################################

#' Rule-of-Thumb Bandwidth Selection
#'
#' Calculates a rule-of-thumb bandwidth for local polynomial regression based on
#' the input matrix `x` and an optional scaling factor `sc`.
#'
#' @param x Numeric matrix of predictor variables.
#' @param sc Numeric, an optional scaling factor for the bandwidth (default is 1).
#'
#' @return Numeric, the calculated bandwidth.
#'
#' @noRd
choose.h = function(x, sc = 1)
{
  n = NROW(x); p = NCOL(x);
  (4 / (2 * p + 1) / n)^(1/(p + 4)) * median(apply(x, 2, sd)) * sc;
}

######################################################################
## THE END
######################################################################
