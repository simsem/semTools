# semTools 0.5-7 (in development)

## New Features:

- `auxiliary()` has 2 new arguments
    - `envir=` is passed to `do.call()`, which prevents problems not finding the `lavaan` package when `semTools` is not loaded via `library()`. Thanks to Brian Keller.
    - `return.syntax=TRUE` will return model syntax for the saturated-correlates parameters, which can be added to the target-model syntax.  This enables adding saturated correlates easily to a `blavaan` model.
- New `fmi(method="cor")` option, to return the fraction missing information for correlation estimates.
- Not necessarily a bug fix, but `MASS::mvrnorm()` has been replaced with `mnormt::rmnorm()` throughout the package.  The latter generates data such that a sample of (e.g.) $N=11$ will have the same first 10 observations as when sampling only $N=10$. 
- New `goricaSEM()` function provided as a wrapper for the `restriktor` package.


## Deprecated Features:

- `efaUnrotate()` deprecated, along with rotation functions `orthRotate()`, `oblqRotate()`, and `funRotate()`.  Users can now directly use `lavaan::efa()` and its `rotation=` argument. Exploratory SEM (ESEM) is also supported by `lavaan` using special `model.syntax()` operators, as shown in its [Issue #112](https://github.com/yrosseel/lavaan/issues/112).
- `runMI()` and all related functionality are deprecated, now returning `OLDlavaan.mi-class` objects.  All `lavaan.mi` functions, class, and methods are now available in the `lavaan.mi` package.  Any `semTools` functions that can operate on either a `lavaan-class` or `lavaan.mi-class` objects (e.g., the `compRelSEM()` function) will still work with the new package's output.


## Bug Fixes:

- Formulas implementing Bollen et al.'s (2012, 2014) BIC extensions have been corrected in the documentation and source code.
- `emmeans` functionality made more robust:
    - to different ordering of names in a specified interaction term
    - to accommodate models fitted to incomplete data using `missing="FIML"` estimation
- `compRelSEM()` no longer gives an error for a `higher=` order factor when indicators are categorical.  A bug for single-group `lavaan.mi` objects was also fixed.
- Fixed an old bug which prevented `clipboard()` from working, particularly on a Mac OS: https://github.com/simsem/semTools/issues/56
- Fixed a bug in the `update()` method for a `measEq.syntax-class` object, preventing an error that occurred when changing values/labels for a parameter across multiple groups.
- Simple intercepts are now always returned by `probe2WayMC()`, `probe2WayRC()`, `probe3WayMC()`, and `probe3WayRC()`, even when the y-intercept is fixed in the fitted model.  This is not actually a bug, but the old default behavior led to `plotProbe()` incorrectly plotting all regression lines through the origin (see https://github.com/simsem/semTools/issues/125 and https://github.com/simsem/semTools/issues/127).
- Fixed a bug in `monteCarloCI()`. See issue [#142](https://github.com/simsem/semTools/issues/142)
- Fixed a bug with `nullRMSEA()` requiring `data=` to be in the global environment; see issue [#133](https://github.com/simsem/semTools/issues/133). Thanks to Brian Keller.
- `lrv2ord()` now works even with a single-variable model.
- `parcelAllocation()` now allows for higher-order factors (see [here](https://groups.google.com/g/lavaan/c/sGb2_9EaeWE/m/xTpxNzGvAAAJ)).


## Known bugs:

- `compareFit()` not currently working for `lavaan.mi-class` objects.



# semTools 0.5-6 (on CRAN 10 May 2022)

## New Features:

- `reliability()` and `reliabilityL2()` deprecated.
    - Average variance extracted now has a dedicated `AVE()` function. 
    - A new `compRelSEM()` function provides more comprehensive options for estimating composite reliability coefficients (various alphas and omegas) than the `reliability()` function.  `compRelSEM()` also estimates reliability for a composite representing a higher-order factor (using the `higher=` argument to name such factors) and implements [Lai's (2021)](https://doi.org/10.1037/met0000287) recently proposed composite-reliability coefficients for different types of multilevel constructs (named with the `config=` and `shared=` arguments).  See https://github.com/simsem/semTools/issues/106 for a conversation with Lai about the implementation, as well as links to documentation of tests.
- `monteCarloCI()` now works for `lavaan.mi` objects, which is much easier than combining multiple imputation with bootstrapping.
- `monteCarloCI()` now optionally returns CIs for the standardized solution, only for `lavaan` objects (not `lavaan.mi`).
- `htmt()` gains an argument `htmt2=TRUE` to use the geometric mean (default) rather than the arithmetic mean (which assumes tau-equivalence).
- `moreFitIndices()` includes 3 new extended BICs discussed by [Bollen et al. (2014)](https://doi.org/10.1177/0049124112452393) 
- `miPowerFit()` can now pass further arguments to `lavaan::modificationIndices()` via `...`
- Added a `tag=` argument to customize character used to flag preferred models in `summary()` method for `FitDiff` class (i.e., `compareFit()` output).  Set `tag=NA` or `NULL` to omit tags.
- `plausibleValues()` applied to `blavaan` objects now passes (via `...`) arguments to `blavaan::blavPredict()`, just as it already passed arguments to `lavPredict()` for objects of `lavaan` or `lavaan.mi` classes.  Thus, setting `type="ypred"` will return plausible values of latent item-responses underlying ordinal observations (analogous to M*plus* option `SAVE: LRESPONSES`).

## Bug Fixes:

- `probe2/3WayM/RC()`: When the moderator was listed first in `nameX=` **and** some (not all) slopes were labeled in the user's model syntax, an error was caused by creating incorrect labels.  Fixes https://github.com/simsem/semTools/issues/103
    - Also, `probe2/3WayM/RC()` help pages have a new **Reference** entry to a tutorial paper written about how to use them.


# semTools 0.5-5 (on CRAN 8 July 2021)

## New Features:

- A new function `lrv2ord()` provides univariate and bivariate population parameters for multivariate-normal data that have been discretized but treated as continuous.
- `reliability()` has new `what=`, `omit.factors=`, and `omit.indicators=` arguments to 
customize composite-reliability estimates.
- `SSpower()` now allows multigroup population models to be specified with `lavaan` syntax via the `popModel=` argument.
- `SSpower()` and `findRMSEApower()` functionality is now available with a graphical user interface:
    - Shiny app available at https://sjak.shinyapps.io/power4SEM/ 
    - tutorial published open-access at https://doi.org/10.3758/s13428-020-01479-0
- `monteCarloMed()` has been replaced by `monteCarloCI()` because the method generalizes beyond evaluating mediation (indirect effects), and can be applied to any function of model parameters. `monteCarloCI()` now accepts a set of functions to sample simultaneously (e.g., both the indirect and total effects), rather than drawing parameters again for each function of parameters. More conveniently, if user-defined parameters are already specified when fitting the `lavaan` model, then `monteCarloCI()` obtains all necessary information from the `lavaan` object.

## Bug Fixes:

- `permuteMeasEq()` returned an error when `parallelType = "multicore"`.
- `reliability()` could use the wrong thresholds when a sequence of ordinal indicator names overlapped (e.g., `x1` and `x10`); now fixed.
- `probe*Way*C()` functions for probing latent interactions formerly returned an error when users specified custom labels for the relevant regression parameters in `lavaan` syntax. And in multigroup models, only group-1 *SE*s were used regardless of which `group=` was specified. Both issues resolved.
- The `show()` method for `compareFit()` now behaves as expected: nothing is printed when assigning `compareFit()` output to an object, but something is printed when there is no assignment. The `summary()` method prints either way and invisibly returns the object itself.


# semTools 0.5-4 (on CRAN 11 January 2021)

## New Features:

- `summary()` method for `lavaan.mi` objects now passes arguments to `lavTestLRT.mi()` via `...`.

## Bug Fixes:

- `reliability()` incorrectly used standardized thresholds rather than the model's parameter estimates (i.e., on the latent-response scales) to calculate omega for categorical indicators.
- `net()` returned an error for models with categorical indicators, now fixed.
- The wrong *df* were used to find confidence limits for RMSEA calculated from the scaled $\chi^2$ statistic in `lavaan.mi` objects.
- `emmeans` methods appropriately deal with incomplete data.
- `probe2WayMC()` and `probe2WayRC()` treated the wrong variable as moderator for point or *SE* estimates, depending on the order in `nameX=`. 


# semTools 0.5-3 (on CRAN 27 May 2020)

## New Functions:

- New `discriminantValidity()` function added.
- Support functions for the `emmeans` package (see `?lavaan2emmeans` for examples)

## New Features:

- `class?lavaan.mi` methods and functions can optionally specify particular imputation numbers to omit, in addition to the general omission criteria in `omit.imps=`. Simply include specific imputation numbers in the `omit.imps=` vector.
- Latent-interaction functions (e.g., `probe2WayMC()`) now work for `lavaan.mi` objects. The `?indProd` documentation now includes an example adding product indicators to multiple imputed data sets.
- The `update()` method for `class?measEq.syntax` has a new argument `change.syntax`. Users can pass lavaan syntax specifying an existing model parameter in order to change the labels or the fixed/free values. This provies some flexibility not found in the `measEq.syntax()` function itself (e.g., releasing an equality constraint in only one of >2 groups, whereas `group.partial=` can only release constraints across all groups).
    - https://github.com/simsem/semTools/issues/60
- The `as.character()` method for `class?measEq.syntax` now accepts the argument `package = "mplus"`, which prints the syntax as an M*plus* MODEL command. The `as.character()` method also has 2 new arguments:
    - `groups.as.blocks=TRUE` optionally prints multigroup syntax in "block" format, which enables users to hack `measEq.syntax()` to specify multilevel CFA models with invariance constraints.
    - `params=` allows users to select specific types of parameters to print, making it easier to check particular aspects of the model specification (e.g., `params = c("loadings","lv.variances")` in metric invariance models).
- `net()` now accepts models fitted to categorical outcomes.
- `reliability()` includes 2 new arguments:
    - default `dropSingle = TRUE` is consistent with old behavior.
    - "total" column no longer returned by default, or ever for 1-factor models. Users can request `return.total = TRUE` when multiple factors are multiple dimensions of a single scale composite.
- New small-*N* corrections added to `chisqSmallN()` function. Also accepts `lavaan.mi` objects now.
- `efaUnrotate()` now accepts summary statistics when `data=NULL`.

## Bug Fixes:

- `reliability()` and `maximalRelia()` returned an error with categorical single-group models
    - https://groups.google.com/d/msg/lavaan/rPVEHUQjqVQ/SQaMrgn-AQAJ
- `reliability()` only ignored higher-order factors without any observed indicators, and returned an error when first-order factors had categorical indicators.  Both issues have been resolved:
    - https://github.com/simsem/semTools/issues/65
- `fitMeasures()` for `lavaan.mi` sometimes returned an error 


# semTools 0.5-2 (on CRAN 30 August 2019)

- Requires `lavaan` version 0.6-5
- Minor bug fixes

- Addition of the `plausibleValues()` function to extract plausible values (i.e., multiple imputations) of factor scores from objects of class `lavaan`, `lavaan.mi`, or `blavaan`
- Full support for `lavaan.mi` models fitted to multiply imputed multilevel data; resolved issue ([#39](https://github.com/simsem/semTools/issues/39)).
- Full support for `fixed.x=TRUE` and `conditional.x=TRUE` in `lavaan.mi` models, including `std.nox` solutions in `summary()`, `modindices.mi()`, and `lavTestScore.mi()`.
- Added the `omit.imps=` argument to all `lavaan.mi`-related functions, optionally excluding solutions that did not converge, failed to estimate *SE*s, or contained NPD matrices (Heywood cases are a special case).  Only the first 2 are excluded by default.
- `reliability()`, `reliabilityL2()`, and `maximalRelia()` now accept `lavaan.mi` objects
- Added (S)EPCs to `lavTestScore.mi()` output when `epc=TRUE`
- Added (A)RIV/FMI to all pooled tests when available for `lavaan.mi` objects, to quantify additional uncertaintly in the test statistic due to missing data.
- Allow multigroup models in `plotProbe()` and related latent-interaction functions.
- `standardizeMx()` was deprecated in previous versions, now removed.



# semTools 0.5-1 (on CRAN 25 September 2018)

- Requires `lavaan` version 0.6-3
- Minor bug fixes
- The formerly deprecated `lisrel2lavaan()` function has been removed from `semTools`

- `compareFit()` now accepts `lavaan.mi` objects returned by `runMI()`
- For `lavaan.mi` objects returned by `runMI()`, the `anova()` method has been updated to behave more like lavaan's `anova()` method:
    - more than 2 nested models can be compared
    - fit indices are no longer an option, and must be requested using the `fitMeasures()` method
- Given the previous addition of score-test functions `modIndices.mi()` and `lavTestScore.mi()` for `lavaan.mi` objects (parallel to `modIndices()` and `lavTestScore()` for `lavaan` objects), the remaining "trilogy" of tests in lavaan (`lavTestLRT()` and `lavTestWald()`) now have parallel functions for `lavaan.mi` objects: `lavTestLRT.mi()` and `lavTestWald.mi()`
    - `lavTestWald.mi()` implements what was formerly available in using `anova(..., test = "D1")`
    - `lavTestLRT.mi()` implements what was formerly available in using `anova(..., test = "D3")`
    - `lavTestLRT.mi()` cannot compare more than 2 nested models. The `anova()` method internally calls the `compareFit()` function to compare multiple `lavaan.mi` objects.
    - For all 3 tests (score, Wald, and LRT), the "D2" pooling method is an option, and there is a newly public function `calculate.D2()` that can be used to pool any set of Wald chi-squared or *z* statistics.
- The `runMI()` function can now be applied to multilevel SEMs that can be fitted with `lavaan()`
    - Known issue ([#39](https://github.com/simsem/semTools/issues/39)): `fitted()` and `resid()` methods will not yet work in models with both multiple levels and multiple groups. This will be resolved in a future version.
- The 3 functions `measurementInvariance()`, `measurementInvarianceCat()`, and `longInvariance()` have been deprecated, redirecting users to the new `measEq.syntax()` function. It is much more general, capable of combining features of all 3 deprecated functions without their restrictions.
    - The function's primary purpose is writing syntax that users can read and edit.
    - The fitting of a model is optional, and fitting multiple models is not (yet) automated. See the `?measEq.syntax` help-page examples for how to fit and compare several levels of invariance.
    - Find many more details [posted on the lavaan forum](https://groups.google.com/d/msg/lavaan/oKwP0_6-i1g/i0jGCU-LBwAJ) (the Google group).



# semTools 0.5-0 (on CRAN 1 July 2018)

- Requires `lavaan` version 0.6-1
- Minor bugs fixed
- Minor convenience features added

- Redesigned `runMI()` function, which no longer produces an object of class `lavaanStar` (that object class is no longer supported), which inherited from class `lavaan`.  It now produces an object of class `lavaan.mi`, which inherits from lavaan's new `lavaanList` class (see the `?lavaanList` help page for details). The reasons to redesign `runMI()` include:
    - The user was required to choose among available methods for pooling the chi-squared test statistic when fitting the model, and the baseline model was also fit so that incremental fit indices (e.g., CFI and TLI) could be calculated.  This lead to more frequent convergence problems.
    - The `lavaanStar` class could inadvertently mislead users into thinking that certain results were available from multiple imputations that were not.  For example, the `modindices()` function would return modification indices for the first imputation, but those were not appropriately pooled statistics.
- The new `runMI()` no longer includes the `chi=` argument, because those options have been moved to an `anova()` method written for `lavaan.mi` objects.  Additional methods have been written: see the `class?lavaan.mi` help page for a list of methods and details about their use.  Additionally, appropriately pooled modification indices and (S)EPCs are now available for multiple imputations (via `modindices.mi()`), as well as a general score test via `lavTestScore.mi()`.
- The `parcelAllocation()` has also been redesigned with new arguments to improve its flexibility and ease of use.  Users are now required to provide lavaan syntax not only for the parcel-level model, but also for the item-level model.  This allows parcel allocation to be automatically detected, no longer requiring users to provide a list of item-to-parcel assignments (see new examples on the help page).
- The `OpenMx` enhancements in previous versions of `semTools` are obsolete now that `OpenMx` provides fit indices and standardized paths, so they have been removed.  However, `standardizeMx()` is still available (temporarily deprecated) to provide standardized mean-structure parameters, until the `OpenMx` maintainers add that feature to `OpenMx::mxStandardizeRAMpaths()`.
- `measurementInvariance()` and `measurementInvarianceCat()` now require users to name all arguments passed to `cfa()`, including the first argument: `model=`.  This is to prevent errors that occurred when some previous users had passed arguments in a different order than expected, which should not be a limitation.


# semTools versions < 0.5

Find our remaining version history on GitHub:

https://github.com/simsem/semTools/wiki/Version-History
