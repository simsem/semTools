# semTools 0.5-2 (submitted to CRAN on 29 August 2019)

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
- The OpenMx enhancements in previous versions of semTools are obsolete now that OpenMx provides fit indices and standardized paths, so they have been removed.  However, `standardizeMx()` is still available (temporarily deprecated) to provide standardized mean-structure parameters, until the OpenMx maintainers add that feature to `OpenMx::mxStandardizeRAMpaths()`.
- `measurementInvariance()` and `measurementInvarianceCat()` now require users to name all arguments passed to `cfa()`, including the first argument: `model=`.  This is to prevent errors that occurred when some previous users had passed arguments in a different order than expected, which should not be a limitation.


# semTools versions < 0.5

Find our remaining version history on GitHub:

https://github.com/simsem/semTools/wiki/Version-History
