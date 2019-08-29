# `semTools`
Useful tools for structural equation modeling.

This is an R package whose primary purpose is to extend the functionality of the R package `lavaan`. There are several *suites* of tools in the package, which correspond to the same theme.  To browse these suites, open the help page at the Console:

```
?semTools::`semTools-package`
```

Additional tools are available to do not require users to rely on any R packages for SEM (e.g., `lavaan`, `OpenMx`, or `sem`), as long as their other software provides the information they need.  Examples:

- `monteCarloMed()` to calculate Monte Carlo confidence intervals for functions of parameters, such as indirect effects in mediation models
- `calculate.D2()` to pool *z* or chi-squared statistics across multiple imputations of missing data
- `indProd()` for creating product indicators of latent interactions
- `SSpower()` provides analytically derived power estimates for SEMs
- `tukeySEM()` for Tukey's WSD post-hoc test of mean-differences under unequal variance and sample size
- `bsBootMiss()` to transform incomplete data to be consistent with the null-hypothesized model, appropriate for model-based (a.k.a. "Bollen--Stine") boostrapping

All users of R (or SEM) are invited to submit functions or ideas for functions by contacting the maintainer, Terrence Jorgensen (TJorgensen314 at gmail dot com). Contributors are encouraged to use **Roxygen** comments to document their contributed code, which is consistent with the rest of `semTools`. Read the vignette from the `roxygen2` package for details: 

```
vignette("rd", package = "roxygen2")
```
