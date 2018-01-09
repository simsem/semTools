
dat <- data.frame(lapply(datCat[1:8], as.numeric))
## no problems with continuous data
mod2con <- efaUnrotate(dat, nf = 2)
summary(mod2con, std = TRUE)
oblqRotate(mod2con)
## try categorical
mod2cat <- efaUnrotate(dat, nf = 2, ordered = paste0("u", 1:8))
summary(mod2cat)
## step through efaUnrotate()
data <- dat
nf <- 2
varList <- colnames(data)
isOrdered <- TRUE
args <- list(ordered = paste0("u", 1:8))



args$data <- data
lavaancfa <- function(...) { lavaan::cfa(...)}
nvar <- length(varList)
facnames <- paste0("factor", 1:nf)
loading <- outer(1:nvar, 1:nf, function(x, y) paste0("load", x, "_", y))
syntax <- ""
for (i in 1:nf) {
variablesyntax <- paste(paste0(loading[,i], "*", varList), collapse = " + ")
factorsyntax <- paste0(facnames[i], " =~ NA*", varList[1], " + ", variablesyntax, "\n")
syntax <- paste(syntax, factorsyntax)
}
syntax <- paste(syntax, paste(paste0(facnames, " ~~ 1*", facnames),
collapse = "\n"), "\n")
syntax
cat(syntax)
if (nf > 1) {
covsyntax <- outer(facnames, facnames,
function(x, y) paste0(x, " ~~ 0*", y, "\n"))[lower.tri(diag(nf), diag = FALSE)]
syntax <- paste(syntax, paste(covsyntax, collapse = " "))
for (i in 2:nf) {
for (j in 1:(i - 1)) {
loadconstraint <- paste(paste0(loading[,i], "*", loading[,j]), collapse=" + ")
syntax <- paste(syntax, paste0("0 == ", loadconstraint), "\n")
}
syntax
loadconstraint
cat(syntax)
args$model <- syntax
do.call(lavaancfa, args)
summary(do.call(lavaancfa, args))
parTable(mod2cat)
syntax <- ""
for (i in 1:nf) {
variablesyntax <- paste(paste0(loading[,i], "*", varList), collapse = " + ")
factorsyntax <- paste0(facnames[i], " =~ NA*", varList[1], " + ", variablesyntax, "\n")
syntax <- paste(syntax, factorsyntax)
}
syntax <- paste(syntax, paste(paste0(facnames, " ~~ 1*", facnames),
collapse = "\n"), "\n")
if (!isOrdered) {
syntax <- paste(syntax, paste(paste0(varList, " ~ 1"), collapse = "\n"), "\n")
}
args$model <- syntax
do.call(lavaancfa, args)
summary(do.call(lavaancfa, args))
covsyntax <- outer(facnames, facnames,
function(x, y) paste0(x, " ~~ 0*", y, "\n"))[lower.tri(diag(nf), diag = FALSE)]
syntax <- paste(syntax, paste(covsyntax, collapse = " "))
args$model <- syntax
summary(do.call(lavaancfa, args))
summary(cfa('f1 =~ u1 + u2 + u3 + u4 ; f2 =~ u5 + u6 + u7 + u8', datCat))
loading
nf
nvar
varList
loading <- outer(1:nvar, 1:nf, function(x, y) paste0("load", x, "_", y))
loading
syntax <- ""
for(i in 1:nf) {
variablesyntax <- paste(paste0(loading[,i], "*", varList), collapse=" + ")
factorsyntax <- paste0(facnames[i], " =~ NA*", varList[1], " + ", variablesyntax, "\n")
syntax <- paste(syntax, factorsyntax)
}
syntax <- paste(syntax, paste(paste0(facnames, " ~~ 1*", facnames), collapse="\n"), "\n")
cat(syntax)
isOrdered
covsyntax <- outer(facnames, facnames, function(x, y) paste0(x, " ~~ 0*", y, "\n"))[lower.tri(diag(nf), diag=FALSE)]
covsyntax
syntax <- paste(syntax, paste(covsyntax, collapse = " "))
cat(syntax)
for(i in 2:nf) {
for(j in 1:(i - 1)) {
loadconstraint <- paste(paste0(loading[,i], "*", loading[,j]), collapse=" + ")
syntax <- paste(syntax, paste0("0 == ", loadconstraint), "\n")
}
cat(syntax)
args <- list(...)
args$model <- syntax
args$data <- data
args
do.call(lavaancfa, args)
summary(do.call(lavaancfa, args))
savehistory("~/svn/semTools/FuturePlans/testEFAcat.R")
