# TensorSIR Package

Welcome to the `TensorSIR` R package, implementing the method proposed in “Parsimonious Tensor Dimension Reduction” by Xing, et al. (2023). This package provides tools for tensor dimensionality reduction and sufficient dimension reduction (SIR) in high-dimensional data analysis.

## Installation
The package is best suited for use on Windows, as it contains binary files that were compiled for Windows.
Clone the package from GitHub:
```{bash}
git clone git@github.com:BioAlgs/TensorSIR.git
```
Then, use devtools to load the package into your R session:
```{r}
devtools::load_all(".")
```
## Usage
Below are some examples demonstrating how to use the `TensorSIR` package. These examples cover basic usage, comparing new and old methodologies, testing hypotheses, and performing cross-validation.
```{r}
# Clean the workspace
rm(list = ls())

# Define matrices
a1 = c(1, 1, 0, 0, 0)
a2 = c(0, 0, 1, 1)
amat = cbind(a1, a2)
b1 = c(1, 1, 1, rep(0, 6))
b2 = c(rep(0, 6), 1, 1, 1)
bmat = cbind(b1, b2)

# Generate data
x = matrix(rnorm(400 * 45), ncol = 45)
boa = kronecker(bmat, amat)
y = (x %*% boa[, 1]) / (2 + (3 + x %*% boa[, 4])^2) + rnorm(400, sd = 0.5)

# Apply TensorSIR
p = c(5, 9)
d = c(1, 1)
nslice = 10
ans1 = foldedSIR(y, x, p, nslice, d)
ans2 = tensorSIR(y, x, p, nslice, 4)

# Compare direction vectors
dist.colspace(boa, ans1$boa)
dist.colspace(boa, ans2$direction)
```
## Comparing methodologies
This example compares the outcomes of the old and new implementations of the foldedSIR function.
```{r}
# Clean the workspace
rm(list = ls())

# [Similar setup as the basic example]

# Old vs. New foldedSIR
ans0 = foldedSIR.old(y, x, p, nslice, c(1, 1))
ans1 = foldedSIR(y, x, p, nslice, c(1, 1))

# Compare column spaces
dist.colspace(ans0$a, ans1$a)
dist.colspace(ans0$b, ans1$b)

boa0 = kronecker(ans0$b, ans0$a)
boa1 = kronecker(ans1$b, ans1$a)
dist.colspace(boa0, boa1)
```

## Testing Hypotheses
Demonstrates how to test hypotheses with the package.
```{r}
# Clean the workspace
rm(list = ls())

# [Data generation steps]

# Perform hypothesis testing
ans = testH0(y, x, nslice)
print(ans$p.value)
wchisq.tail(ans$stat, ans$z)
```

## Cross-Validation
Example showing how to perform cross-validation.
```{r}
# Clean the workspace
rm(list = ls())

# [Data generation steps]

# Cross-validation
ans = CV4B(y, x, boa[, c(1, 4)], 0.5)
```
## Conclusion
We hope this package assists in your research and analysis. For any issues or contributions, please open an issue or pull request on GitHub.

