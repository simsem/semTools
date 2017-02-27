\name{auxiliary}
\alias{auxiliary}
\alias{cfa.auxiliary}
\alias{sem.auxiliary}
\alias{growth.auxiliary}
\alias{lavaan.auxiliary}
\title{
	Implement Saturated Correlates with FIML
}
\description{
	Automatically add auxiliary variables to a lavaan model when using full information maximum likelihood (FIML) to handle missing data
}
\usage{
auxiliary(model, data, aux, fun, ...)
lavaan.auxiliary(model, data, aux, ...)
cfa.auxiliary(model, data, aux, ...)
sem.auxiliary(model, data, aux, ...)
growth.auxiliary(model, data, aux, ...)
}
\arguments{
  \item{model}{
	The analysis model can be specified with 1 of 2 objects:
   \enumerate{
     \item  lavaan \code{\link[lavaan]{model.syntax}} specifying a hypothesized
            model \emph{without} mention of auxiliary variables in \code{aux}
     \item  a parameter table, as returned by \code{\link[lavaan]{parTable}},
            specifying the target model \emph{without} auxiliary variables.
            This option requires these columns (and silently ignores all others):
            \code{c("lhs","op","rhs","user","group","free","label","plabel","start")}
   }
}
  \item{data}{
	\code{data.frame} that includes auxiliary variables as well as any observed variables in the \code{model}
}
  \item{aux}{
	\code{character}. Names of auxiliary variables to add to \code{model}
  }
  \item{fun}{
	\code{character}. Name of a specific lavaan function used to fit \code{model} to \code{data} (i.e., \code{"lavaan"}, \code{"cfa"}, \code{"sem"}, or \code{"growth"}). Only required for \code{auxiliary}.
  }
  \item{\dots}{
	additional arguments to pass to \code{\link[lavaan]{lavaan}}. See also \code{\link[lavaan]{lavOptions}}
  }
}
\details{
These functions are wrappers around the corresponding lavaan functions.
You can use them the same way you use \code{\link[lavaan]{lavaan}}, but you
\emph{must} pass your full \code{data.frame} to the \code{data} argument.
Because the saturated-correlates approaches (Enders, 2008) treates exogenous
variables as random, \code{fixed.x} must be set to \code{FALSE}. Because FIML
requires continuous data (although nonnormality corrections can still be requested),
no variables in the model nor auxiliary variables specified in \code{aux}
can be declared as \code{ordered}.
}
\value{
	a fitted \code{\linkS4class{lavaan}} object.  Additional
   information is stored as a \code{list} in the \code{@external} slot:
   \itemize{
     \item \code{$baseline.model}. a fitted \code{\linkS4class{lavaan}}
           object. Results of fitting an appropriate independence model for
           the calculation of incremental fit indices (e.g., CFI, TLI) in
           which the auxiliary variables remain saturated, so only the target
           variables are constrained to be orthogonal. See \strong{Examples} for how
           to send this baseline model to \code{\link[lavaan]{fitMeasures}}.
     \item \code{$aux}. The character vector of auxiliary variable names.
   }
}
\references{
Enders, C. K. (2008). A note on the use of missing auxiliary variables in full information maximum likelihood-based structural equation models. \emph{Structural Equation Modeling, 15}(3), 434-448. doi:10.1080/10705510802154307
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\examples{
HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
                                      "ageyr","agemo","school")]
set.seed(12345)
HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
age <- HSMiss$ageyr + HSMiss$agemo/12
HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)

## specify CFA model from lavaan's ?cfa help page
HS.model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
'

## works just like cfa(), but with an extra "aux" argument
fitaux1 <- cfa.auxiliary(HS.model, data = HSMiss, aux = c("ageyr","agemo"),
                         missing = "fiml", estimator = "mlr")

## example with multiple groups
fitaux2 <- cfa.auxiliary(HS.model, data = HSMiss, aux = c("ageyr","agemo"),
                         group = "school", group.equal = "loadings")

}
