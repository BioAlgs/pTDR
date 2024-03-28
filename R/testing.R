######################################################################
## test the number of directions
## 02-06-2015
## Peng Zeng
######################################################################
#' @useDynLib testing-cfuns 
#' @importFrom dr dr.slices
#' @export
dr::dr.slices

######################################################################
## H0: Mmat = 0
## test statistic: tr(Mmat) 
######################################################################

#' Test Hypothesis for Independence in SIR Model
#'
#' Performs a hypothesis test for independence between response `y` and predictors `x`
#' within the SIR model framework, using the number of slices specified by `nslice`.
#'
#' @param y Numeric vector of response variables.
#' @param x Numeric matrix of predictor variables.
#' @param nslice Integer specifying the number of slices for analysis.
#'
#' @return A list containing various statistics including the test statistic (`stat`), 
#'   p-value (`p.value`), eigenvalues (`z`), degrees of freedom (`df`), and scale factor (`scale`).
#'
#' @examples
#' y <- rnorm(100)
#' x <- matrix(rnorm(200), ncol = 2)
#' nslice <- 10
#' result <- testH0(y, x, nslice)
#' print(result)
#'
#' @export
testH0 = function(y, x, nslice)
{
  n = length(y);

  temp = SIR.slice(y, x, nslice)
  stat = n * sum(diag(temp$Mmat));

  info = dr.slices(y, nslice);
  psimat = C.calpsimat(x, info);
  z = eigen(psimat, only.values = TRUE)$values / n; 

  zsum = mean(diag(psimat));
  z2sum = mean(psimat * psimat);

  ## approximate \sum z_k \chi_1^2 by g \chi_f^2
  ## evaluate tail probability
  ## where g = sum(z^2) / sum(z), f = [sum(z)]^2 / sum(z^2)

  chisq.scale = z2sum / zsum;
  chisq.df = zsum * zsum / z2sum;
  chisq.pval = pchisq(stat / chisq.scale, df = chisq.df, lower.tail = FALSE);

  list(sliceinfo = info, psimat = psimat, 
       stat = stat, p.value = chisq.pval, z = z,
       df = chisq.df, scale = chisq.scale);
}

######################################################################
## evaluate p-value for weighted chi-squared distribution
######################################################################

#' Evaluate P-Value for Weighted Chi-Squared Distribution
#'
#' Calculates the tail probability of a weighted chi-squared distribution,
#' given a test statistic, eigenvalues vector `zvec`, and the number of Monte Carlo simulations.
#'
#' @param stat Numeric, the test statistic for which to evaluate the tail probability.
#' @param zvec Numeric vector of eigenvalues used in the weighting.
#' @param n.mc Integer, the number of Monte Carlo simulations to perform (default = 10000).
#' @param eps Numeric, a small tolerance value to filter `zvec` elements (default = 1e-7).
#'
#' @return Numeric, the estimated tail probability.
#'
#' @export
wchisq.tail = function(stat, zvec, n.mc = 10000, eps = 1e-7)
{
  index = (zvec > eps);
  m = sum(index);
  z = zvec[index];
  ans = numeric(m);
  for(i in 1:n.mc)
  {
    ans[i] = sum(rchisq(m, df = 1) * z)
  }

  mean(ans > stat);
}

######################################################################
## calculate matrix Mmat, psi(xi, xj)
######################################################################

#' Calculate Matrix Psi for SIR Model
#'
#' Computes the psi matrix used in SIR model analysis, based on predictors `x`
#' and slice information.
#'
#' @param x Numeric matrix of predictor variables.
#' @param sliceinfo List containing slice sizes and indicators.
#'
#' @return Numeric matrix, the computed psi matrix.
#'
#' @noRd
calpsimat = function(x, sliceinfo)
{
  n = NROW(x);
  psimat = matrix(0, nrow = n, ncol = n);

  prop = sliceinfo$slice.sizes / n;

  for(i in 1:n)
  for(j in i:n)
  {
    hi = sliceinfo$slice.indicator[i];
    hj = sliceinfo$slice.indicator[j];

    psi = -sum(x[i, ] * x[j, ]);

    if(hi == hj)
    {
      psi = psi * (1.0 - 1.0 / prop[hi]);
    }

    psimat[i, j] = psi;
    psimat[j, i] = psi;
  }

  psimat;
}

#' Interface to C Function for Calculating Psi Matrix
#'
#' Provides an interface to a C function for calculating the psi matrix, optimizing 
#' the computation compared to pure R implementation.
#'
#' @param x Numeric matrix of predictor variables.
#' @param sliceinfo List containing slice sizes and indicators.
#'
#' @return Numeric matrix, the psi matrix calculated by the C function.
#'
#' @noRd
C.calpsimat = function(x, sliceinfo)
{
  n = NROW(x); p = NCOL(x);
  psimat = matrix(0, nrow = n, ncol = n);

  prop = sliceinfo$slice.sizes / n;

  ans = .C("calpsimat", as.integer(n), as.integer(p), as.double(x), 
           as.integer(sliceinfo$slice.indicator - 1), as.double(prop), 
           psi = double(n*n),
           PACKAGE = "testing-cfuns");

  matrix(ans$psi, nrow = n, ncol = n);
}

######################################################################
## THE END
######################################################################
