\dontrun{

  library(lavaan)
  library(emmeans)

  #### Moderation Analysis ####

  mean_sd <- function(x) mean(x) + c(-sd(x), 0, sd(x))

  model1 <- '
  # regressions
  Sepal.Length ~ b1 * Sepal.Width + b2 * Petal.Length + b3 * Sepal.Width:Petal.Length


  # define mean parameter label for centered math for use in simple slopes
  Sepal.Width ~ Sepal.Width.mean * 1

  # define variance parameter label for centered math for use in simple slopes
  Sepal.Width ~~ Sepal.Width.var * Sepal.Width

  # simple slopes for condition effect
  SD.below := b2 + b3 * (Sepal.Width.mean - sqrt(Sepal.Width.var))
  mean     := b2 + b3 * (Sepal.Width.mean)
  SD.above := b2 + b3 * (Sepal.Width.mean + sqrt(Sepal.Width.var))
  '

  semFit1 <- sem(model = model1,
                 data = iris)

  ## Compare simple slopes
  # From `emtrends`
  test(
    emtrends(semFit1, ~ Sepal.Width, "Petal.Length",
             lavaan.DV = "Sepal.Length",
             cov.red = mean_sd)
  )

  # From lavaan
  parameterEstimates(semFit1, output = "pretty")[13:15, ]
  # Identical slopes.
  # SEs differ due to lavaan estimating uncertainty of the mean / SD
  # of Sepal.Width, whereas emmeans uses the mean+-SD as a is.


  ## Compare emmeans
  # From lm -> emmeans
  lmFit <- lm(Sepal.Length ~ Sepal.Width * Petal.Length, iris)
  emmeans(lmFit, ~ Sepal.Width, at = list(Sepal.Width = 1:2))

  # From lavaan -> emmeans
  emmeans(semFit1, ~ Sepal.Width, at = list(Sepal.Width = 1:2),
          lavaan.DV = "Sepal.Length")

  #### Latent DV ####

  model2 <- '
  LAT1 =~ Sepal.Length + Sepal.Width

  LAT1 ~ b1 * Petal.Width + 1 * Petal.Length

  Petal.Length ~ Petal.Length.mean * 1

  V1 := 1 * Petal.Length.mean + 1 * b1
  V2 := 1 * Petal.Length.mean + 2 * b1
  '

  semFit2 <- sem(model = model2,
                 data = iris, std.lv = TRUE)

  ## Compare emmeans
  # From emmeans
  test(
    emmeans(semFit2, ~ Petal.Width,
            lavaan.DV = "LAT1",
            at = list(Petal.Width = 1:2))
  )

  # From lavaan
  parameterEstimates(semFit2, output = "pretty")[15:16, ]
  # Identical means.
  # SEs differ due to lavaan estimating uncertainty of the mean
  # of Petal.Length, whereas emmeans uses the mean as is.

  #### Dealing with factors ####

  warpbreaks <- cbind(warpbreaks,
                      model.matrix(~ wool + tension, data = warpbreaks))

  model3 <- "
  # Split for convenience
  breaks ~ 1
  breaks ~ woolB
  breaks ~ tensionM + tensionH
  breaks ~ woolB:tensionM + woolB:tensionH
  "

  semFit3 <- sem(model3, warpbreaks)

  ## Compare contrasts
  # From lm -> emmeans
  lmFit2 <- lm(breaks ~ wool * tension, data = warpbreaks)
  lmEM <- emmeans(lmFit2, ~ tension + wool)
  contrast(lmEM, method = data.frame(L_all = c(-1, .05, 0.5),
                                     M_H   = c(0, 1, -1)), by = "wool")

  # From lavaan -> emmeans
  lavEM <- emmeans(semFit3, ~ tensionM + tensionH + woolB,
                   lavaan.DV = "breaks")
  contrast(lavEM,
           method = data.frame(
             "L_all|A" = c(c(-1, .05, 0.5, 0), rep(0, 4)),
             "M_H|A"   = c(c(0, 1, -1, 0),     rep(0, 4)),
             "L_all|A" = c(rep(0, 4),          c(-1, .05, 0.5, 0)),
             "M_H|A"   = c(rep(0, 4),          c(0, 1, -1, 0))
           ))
}

