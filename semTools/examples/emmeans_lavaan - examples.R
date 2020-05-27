\dontrun{

  library(lavaan)
  library(emmeans)

  #### Moderation Analysis ####

  mean_sd <- function(x) mean(x) + c(-sd(x), 0, sd(x))

  model <- '
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

  semFit <- sem(model = model,
                data = iris)

  ## Compare simple slopes
  # From `emtrends`
  test(
    emtrends(semFit, ~ Sepal.Width, "Petal.Length",
             lavaan.DV = "Sepal.Length",
             cov.red = mean_sd)
  )

  # From lavaan
  parameterEstimates(semFit, output = "pretty")[13:15, ]
  # Identical slopes.
  # SEs differ due to lavaan estimating uncertainty of the mean / SD
  # of Sepal.Width, whereas emmeans uses the mean+-SD as is (fixed).


  #### Latent DV ####

  model <- '
  LAT1 =~ Sepal.Length + Sepal.Width

  LAT1 ~ b1 * Petal.Width + 1 * Petal.Length

  Petal.Length ~ Petal.Length.mean * 1

  V1 := 1 * Petal.Length.mean + 1 * b1
  V2 := 1 * Petal.Length.mean + 2 * b1
  '

  semFit <- sem(model = model,
                data = iris, std.lv = TRUE)

  ## Compare emmeans
  # From emmeans
  test(
    emmeans(semFit, ~ Petal.Width,
            lavaan.DV = "LAT1",
            at = list(Petal.Width = 1:2))
  )

  # From lavaan
  parameterEstimates(semFit, output = "pretty")[15:16, ]
  # Identical means.
  # SEs differ due to lavaan estimating uncertainty of the mean
  # of Petal.Length, whereas emmeans uses the mean as is.

  #### Multi-Variate DV ####

  model <- '
  ind60 =~ x1 + x2 + x3

  # metric invariance
  dem60 =~ y1 + a*y2 + b*y3 + c*y4
  dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # scalar invariance
  y1 + y5 ~ d*1
  y2 + y6 ~ e*1
  y3 + y7 ~ f*1
  y4 + y8 ~ g*1

  # regressions (slopes differ: interaction with time)
  dem60 ~ b1*ind60
  dem65 ~ b2*ind60 + NA*1 + Mean.Diff*1

  # residual correlations
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8

  # conditional mean differences (besides mean(ind60) == 0)
   low := (-1*b2 + Mean.Diff) - (-1*b1) # 1 SD below M
  high := (b2 + Mean.Diff) - b1         # 1 SD above M
'

  semFit <- sem(model, data = PoliticalDemocracy)


  ## Compare contrasts
  # From emmeans
  emmeans(semFit, pairwise ~ rep.meas|ind60,
          lavaan.DV = c("dem60","dem65"),
          at = list(ind60 = c(-1,1)))[[2]]

  # From lavaan
  parameterEstimates(semFit, output = "pretty")[49:50, ]


  #### Multi Group ####

  model <- 'x1 ~ c(int1, int2)*1 + c(b1, b2)*ageyr
  diff_11 := (int2 + b2*11) - (int1 + b1*11)
  diff_13 := (int2 + b2*13) - (int1 + b1*13)
  diff_15 := (int2 + b2*15) - (int1 + b1*15)
'
  semFit <- sem(model, group = "school", data = HolzingerSwineford1939)


  ## Compare contrasts
  # From emmeans (note `nesting = NULL`)
  emmeans(semFit, pairwise ~ school | ageyr, lavaan.DV = "x1",
          at = list(ageyr = c(11, 13, 15)), nesting = NULL)[[2]]

  # From lavaan
  parameterEstimates(semFit, output = "pretty")

  #### Dealing with factors ####

  warpbreaks <- cbind(warpbreaks,
                      model.matrix(~ wool + tension, data = warpbreaks))

  model <- "
  # Split for convenience
  breaks ~ 1
  breaks ~ woolB
  breaks ~ tensionM + tensionH
  breaks ~ woolB:tensionM + woolB:tensionH
  "

  semFit <- sem(model, warpbreaks)

  ## Compare contrasts
  # From lm -> emmeans
  lmFit <- lm(breaks ~ wool * tension, data = warpbreaks)
  lmEM <- emmeans(lmFit, ~ tension + wool)
  contrast(lmEM, method = data.frame(L_all = c(-1, .05, 0.5),
                                     M_H   = c(0, 1, -1)), by = "wool")

  # From lavaan -> emmeans
  lavEM <- emmeans(semFit, ~ tensionM + tensionH + woolB,
                   lavaan.DV = "breaks")
  contrast(lavEM,
           method = list(
             "L_all|A" = c(c(-1, .05, 0.5, 0), rep(0, 4)),
             "M_H  |A" = c(c(0, 1, -1, 0),     rep(0, 4)),
             "L_all|A" = c(rep(0, 4),          c(-1, .05, 0.5, 0)),
             "M_H  |A" = c(rep(0, 4),          c(0, 1, -1, 0))
           ))
}
