######################################################################
## local polynomial smoothing
######################################################################

dyn.load(paste("./inst/libs/locpoly", .Platform$dynlib.ext, sep = ""));

Clocpoly = function(yin, xin, h, xout)
{
  n = length(yin); p = NCOL(xin); m = NROW(xout);
  if(n != NROW(xin))  stop('NROW(xin) != length(yin).\n');
  if(p != NCOL(xout)) stop('NCOL(xin) != NCOL(xout).\n');

  .C('locpoly', as.integer(n), as.integer(p), as.integer(m), 
      as.double(yin), as.double(xin), as.double(h),
      as.double(xout), yhat = double(m))$yhat;
}


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

choose.h = function(x, sc = 1)
{
  n = NROW(x); p = NCOL(x);
  (4 / (2 * p + 1) / n)^(1/(p + 4)) * median(apply(x, 2, sd)) * sc;
}

######################################################################
## THE END
######################################################################
