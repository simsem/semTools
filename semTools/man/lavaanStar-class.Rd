\name{lavaanStar-class}
\docType{class}
\alias{lavaanStar-class}
\alias{inspect,lavaanStar-method}
\alias{summary,lavaanStar-method}

\title{Class For Representing A (Fitted) Latent Variable Model with Additional Elements}
\description{This is the \code{lavaan} class that contains additional information about the fit values from the null model. Some functions are adjusted according to the change.}
\section{Objects from the Class}{
Objects can be created via the \code{\link{auxiliary}} function or \code{\link{runMI}}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{The function call as returned by \code{match.called()}.}
    \item{\code{timing}:}{The elapsed time (user+system) for various parts of 
      the program as a list, including the total time.}
    \item{\code{Options}:}{Named list of options that were provided by
      the user, or filled-in automatically.}
    \item{\code{ParTable}:}{Named list describing the model parameters. Can be coerced to a data.frame. In the documentation, this is called the `parameter table'.}
    \item{\code{Data}:}{Object of internal class \code{"Data"}: information
about the data.}
    \item{\code{SampleStats}:}{Object of internal class \code{"SampleStats"}: sample 
      statistics}
    \item{\code{Model}:}{Object of internal class \code{"Model"}: the 
      internal (matrix) representation of the model}
    \item{\code{Fit}:}{Object of internal class \code{"Fit"}: the 
      results of fitting the model}
	\item{\code{nullfit}:}{The fit-indices information from the null model}
	\item{\code{imputed}:}{The list of information from running multiple imputation. The first element is the convergence rate of the target and null models. The second element is the fraction missing information. The first estimate of FMI (FMI.1) is asymptotic FMI and the second estimate of FMI (FMI.2) is corrected for small numbers of imputation. The third element is the fit values of the target model by the specified chi-squared methods. The fourth element is the fit values of the null model by the specified chi-square methods.}
	\item{\code{auxNames}:}{The list of auxiliary variables in the analysis.}
  }
}
\references{
see \code{\linkS4class{lavaan}}}
\seealso{
\code{\link{auxiliary}}; \code{\link{runMI}}
}
\author{
    Sunthud Pornprasertmanit (Texas Tech University; \email{sunthud.pornprasertmanit@ttu.edu})
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
			  
dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
			  
fit <- cfa(HS.model, data=dat) 
fitaux <- auxiliary(HS.model, aux="z", data=dat, fun="cfa")
}
