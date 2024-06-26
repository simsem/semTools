% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measEq.R
\docType{class}
\name{measEq.syntax-class}
\alias{measEq.syntax-class}
\alias{show,measEq.syntax-method}
\alias{summary,measEq.syntax-method}
\alias{as.character,measEq.syntax-method}
\alias{update,measEq.syntax-method}
\title{Class for Representing a Measurement-Equivalence Model}
\usage{
\S4method{as.character}{measEq.syntax}(x, package = "lavaan",
  params = NULL, single = TRUE, groups.as.blocks = FALSE)

\S4method{show}{measEq.syntax}(object)

\S4method{summary}{measEq.syntax}(object, verbose = TRUE)

\S4method{update}{measEq.syntax}(object, ..., evaluate = TRUE,
  change.syntax = NULL)
}
\arguments{
\item{x, object}{an object of class \code{measEq.syntax}}

\item{package}{\code{character} indicating the package for which the model
syntax should be generated.  Currently, only \code{"lavaan"} and
\code{"mplus"} are supported.}

\item{params}{\code{character} vector indicating which type(s) of parameter
to print syntax for. Must match a type that can be passed to
\code{group.equal} or \code{long.equal}, but \code{"residual.covariances"}
and \code{"lv.covariances"} will be silently ignored. Instead, requesting
\code{"residuals"} or \code{"lv.variances"} will return covariances along
with variances. By default (\code{NULL}), all types are printed.}

\item{single}{\code{logical} indicating whether to concatenate lavaan
\code{\link[lavaan:model.syntax]{lavaan::model.syntax()}} into a single \code{character} string.
Setting \code{FALSE} will return a vector of strings, which may be
convenient (or even necessary to prevent an error) in
models with long variable names, many variables, or many groups.}

\item{groups.as.blocks}{\code{logical} indicating whether to write lavaan
\code{\link[lavaan:model.syntax]{lavaan::model.syntax()}} using vectors of labels and values
for multiple groups (the default: \code{FALSE}), or whether to write
a separate "block" of syntax per group. The block structure could allow
users to apply the generated multigroup syntax (after some editing) to
test invariance across levels in a multilevel SEM (see final example on
\code{\link[=measEq.syntax]{measEq.syntax()}} help page).}

\item{verbose}{\code{logical} indicating whether to print a summary to the
screen (default). If \code{FALSE}, only a pattern matrix is returned.}

\item{...}{Additional arguments to the \code{call}, or arguments with
changed values.}

\item{evaluate}{If \code{TRUE}, evaluate the new \code{call}; otherwise,
return the new \code{call}.}

\item{change.syntax}{\code{\link[lavaan:model.syntax]{lavaan::model.syntax()}} specifying
labels or fixed/free values of parameters in \code{object}.
These provide some flexibility to customize
existing parameters without having to copy/paste the output of
\code{as.character(object)} into an R script. For example,
\code{group.partial} will free a parameter across all groups, but
\code{update} allows users to free the parameter in just one group
while maintaining equality constraints among other groups.}
}
\value{
\item{summary}{\code{signature(object = "measEq.syntax", verbose = TRUE)}:
A \code{character} matrix indicating the pattern of \code{numeric},
\code{ordered}, or latent indicators loading on common factors.
By default (\code{verbose = TRUE}), \code{summary} also prints descriptive
details about the model, including the numbers of indicators and factors,
and which parameters are constrained to equality.}
\item{show}{\code{signature(object = "measEq.syntax")}: Prints a message
about how to use the \code{object} for model fitting. Invisibly
returns the \code{object}.}
\item{update}{\code{signature(object = "measEq.syntax", ..., evaluate = TRUE, change.syntax = NULL)}: Creates a new
\code{object} with updated arguments in \code{...}, or updated
parameter labels or fixed/free specifications in \code{object}.}
\item{as.character}{\code{signature(x = "measEq.syntax", package = "lavaan")}:
Converts the \code{measEq.syntax} object to model syntax that can be
copy/pasted or written to a syntax file to be edited before analysis,
or simply passed to \code{\link[lavaan:lavaan]{lavaan::lavaan()}} to fit the model to
data. Generated M\emph{plus} syntax could also be utilized using the
\pkg{MplusAuthomation} package.}
}
\description{
This class of object stores information used to automatically generate
lavaan model syntax to represent user-specified levels of measurement
equivalence/invariance across groups and/or repeated measures. See
\code{\link[=measEq.syntax]{measEq.syntax()}} for details.
}
\section{Slots}{

\describe{
\item{\code{package}}{\code{character} indicating the software package used to
represent the model. Currently, only \code{"lavaan"} is available, which
uses the LISREL representation (see \code{\link[lavaan:lavOptions]{lavaan::lavOptions()}}).
In the future, \code{"OpenMx"} may become available, using RAM
representation.}

\item{\code{model.type}}{\code{character}. Currently, only "cfa" is available.
Future versions may allow for MIMIC / RFA models, where invariance can be
tested across levels of exogenous variables explicitly included as
predictors of indicators, controlling for their effects on (or correlation
with) the common factors.}

\item{\code{call}}{The function call as returned by \code{match.call()}, with
some arguments updated if necessary for logical consistency.}

\item{\code{meanstructure}}{\code{logical} indicating whether a mean structure is
included in the model.}

\item{\code{numeric}}{\code{character} vector naming \code{numeric} manifest indicators.}

\item{\code{ordered}}{\code{character} vector naming \code{ordered} indicators.}

\item{\code{parameterization}}{\code{character}. See \code{\link[lavaan:lavOptions]{lavaan::lavOptions()}}.}

\item{\code{specify}}{\code{list} of parameter matrices, similar in form to the
output of \code{lavInspect(fit, "free")}. These matrices
are \code{logical}, indicating whether each parameter should be specified
in the model syntax.}

\item{\code{values}}{\code{list} of parameter matrices, similar in form to the
output of \code{lavInspect(fit, "free")}. These matrices
are \code{numeric}, indicating whether each parameter should be freely
estimated (indicated by \code{NA}) or fixed to a particular value.}

\item{\code{labels}}{\code{list} of parameter matrices, similar in form to the
output of \code{lavInspect(fit, "free")}. These matrices
contain \code{character} labels used to constrain parameters to equality.}

\item{\code{constraints}}{\code{character} vector containing additional equality
constraints used to identify the model when \code{ID.fac = "fx"}.}

\item{\code{ngroups}}{\code{integer} indicating the number of groups.}
}}

\examples{
## See ?measEq.syntax help page for examples using lavaan


}
\author{
Terrence D. Jorgensen (University of Amsterdam;
\email{TJorgensen314@gmail.com})
}
