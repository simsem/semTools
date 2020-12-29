### Terrence D. Jorgensen
### Last updated: 30 October 2019
### new DIF-detection tools for measEq suite


findAnchors <- function(object, test = c("LRT","Wald","score"),
                        nAnchors = 2L, # or a proportion
                        method = "woods.2009", # or kopf's iterative routine(s)?
                        criterion, # lowest chi-squared, alpha level, effect size?
                        conParams = c("thresholds","loadings"."intercepts"),   # like param= from permuteMeasEq()
                        exclude = "", # exceptions (variables to test for DIF, not anchors)
                        longFacNames = list(), # as in measEq.syntax
                        longIndNames = list()) {

}



