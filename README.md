# pTDR R Package

Welcome to the `pTDR` R package, implementing the method proposed in “Parsimonious Tensor Dimension Reduction” by Xing, et al. (2023). This package provides tools for tensor dimensionality reduction and sufficient dimension reduction (SIR) in high-dimensional data analysis.

## Requirements

- Windows System
- R >= 4.3.1

## Installation
The package is best suited for use on Windows, as it contains binary files that were compiled for Windows.
Clone the package from GitHub:
```{bash}
git clone https://github.com/BioAlgs/pTDR.git
```
Then, in start R or Rstudio, set the working directory as the downloaded folder ("./pTDR") and use devtools to load the package into your R session:
```{r}
devtools::load_all(".")
```
## Help File

You can access the help document of main functions such as `pTDR`, `testH0` by using `? pTDR` and `? testH0`.

The help file of `pTDR` is here:

**Usage**
```
pTDR(y, x, p, nslice, ndir = 1, maxiter = 200, tol = 1e-08, eps = 1e-08)
```
**Arguments**

`y`	Numeric vector of response variables, indicating the outcome for each observation.

`x`	Numeric matrix of predictor variables, where each row is an observation and columns correspond to variables, possibly representing a flattened tensor.

`p`	Numeric vector indicating the dimensions of the tensor representation for the predictor variables. The product of p should equal the number of columns in x.

`nslice` Integer specifying the number of slices for the SIR analysis.

`ndir` Integer indicating the number of directions to estimate (default is 1).

`maxiter` Maximum number of iterations for the algorithm (default is 200).

`tol`	Tolerance level for convergence in the optimization algorithm (default is 1e-8).

`eps`	A small numerical value to stabilize computations (default is 1e-8).

**Value**

A list containing the following components:

- `p`: The dimensions of the tensor as specified by the input.

- `ndir`: The number of directions estimated.

- `info`: Information from the initial SIR analysis, including slice information.

- `ans`: Results from the tensor-based SIR analysis, varying based on the dimensionality of p.

## Examples

Below are some examples demonstrating how to use the `pTDR` package. These examples cover basic usage, comparing new and old methodologies, testing hypotheses, and performing cross-validation.
```{r}
# Define matrices
a1 = c(1, 1, 0, 0, 0)
a2 = c(0, 0, 0, 1, 1)
amat = cbind(a1, a2)
b1 = c(1, 1, 1, rep(0, 6))
b2 = c(rep(0, 6), 1, 1, 1)
bmat = cbind(b1, b2)

# Generate data
x = matrix(rnorm(400 * 45), ncol = 45)
boa = kronecker(bmat, amat)
y = (x %*% boa[, 1]) / (2 + (3 + x %*% boa[, 4])^2) + rnorm(400, sd = 0.5)

# Apply pTDR
p = c(5, 9)
d = c(1, 1)
nslice = 10
ans1 = foldedSIR(y, x, p, nslice, d)
ans2 = pTDR(y, x, p, nslice, 4)

# Compare direction vectors
dist.colspace(boa, ans1$boa)
dist.colspace(boa, ans2$direction)
```


## Testing Hypotheses
Demonstrates how to test hypotheses with the package.
```{r}
# [Data generation steps]

# Perform hypothesis testing
ans = testH0(y, x, nslice)
print(ans$p.value)
wchisq.tail(ans$stat, ans$z)
```



## Conclusion
We hope this package assists in your research and analysis. For any issues or contributions, please open an issue or pull request on GitHub.

