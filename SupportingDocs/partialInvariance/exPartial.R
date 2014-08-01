library(semTools)

HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

models2 <- measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school", std.lv = TRUE)

# Note that scalar invariance cannot be established

# Build-up strategy

round(partialInvariance(models2, type = "scalar", p.adjust = "holm"), 4)
round(partialInvariance(models2, type = "scalar", fix = "x1", p.adjust = "holm"), 4)
round(partialInvariance(models2, type = "scalar", fix = c("x1", "x2"), p.adjust = "holm"), 4)
round(partialInvariance(models2, type = "scalar", fix = c("x1", "x2", "x9"), p.adjust = "holm"), 4)

result <- partialInvariance(models2, type = "scalar", fix = c("x1", "x2", "x9"), p.adjust="holm", return.fit = TRUE)
summary(result$models$parent)

# Tear-down strategy

round(partialInvariance(models2, type = "scalar", p.adjust = "holm"), 4)
round(partialInvariance(models2, type = "scalar", free="x3", p.adjust = "holm"), 4)
round(partialInvariance(models2, type = "scalar", free=c("x3", "x7"), p.adjust = "holm"), 4)

result <- partialInvariance(models2, type = "scalar", free=c("x3", "x7"), p.adjust="holm", return.fit = TRUE)
summary(result$models$nested)


