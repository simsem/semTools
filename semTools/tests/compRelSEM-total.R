library(semTools)

cfa_data <- lavaan::HolzingerSwineford1939
higher_order <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
  general =~ visual + textual + speed
"
all_items <- ".TOTAL. <~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"

for (group in list(NULL, "school")) {
  fit <- lavaan::cfa(higher_order, data = cfa_data, group = group)
  with_total <- compRelSEM(fit, return.total = TRUE)
  explicit_total <- compRelSEM(fit, W = all_items)

  stopifnot(
    ".TOTAL." %in% names(with_total),
    isTRUE(all.equal(with_total[".TOTAL."], explicit_total)),
    isTRUE(all.equal(compRelSEM(fit, return.total = -1L), explicit_total))
  )
}

# x9 is cross-loaded before the remaining speed indicators.
cross_loading <- "
  visual =~ x1 + x2 + x3 + x9
  speed  =~ x9 + x7 + x8
"
fit <- lavaan::cfa(cross_loading, data = cfa_data)
total <- compRelSEM(fit, return.total = -1L)
explicit_total <- compRelSEM(fit, W = ".TOTAL. <~ x1 + x2 + x3 + x7 + x8 + x9")
stopifnot(isTRUE(all.equal(total, explicit_total)))
