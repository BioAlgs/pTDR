######################################################################
## test the number of directions
## 02-06-2015
## Peng Zeng
######################################################################

dyn.load(paste("./inst/libs/testing-cfuns", .Platform$dynlib.ext, sep = ""));

######################################################################
## H0: Mmat = 0
## test statistic: tr(Mmat) 
######################################################################

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


C.calpsimat = function(x, sliceinfo)
{
  n = NROW(x); p = NCOL(x);
  psimat = matrix(0, nrow = n, ncol = n);

  prop = sliceinfo$slice.sizes / n;

  ans = .C("calpsimat", as.integer(n), as.integer(p), as.double(x), 
           as.integer(sliceinfo$slice.indicator - 1), as.double(prop), 
           psi = double(n*n));

  matrix(ans$psi, nrow = n, ncol = n);
}

######################################################################
## THE END
######################################################################
