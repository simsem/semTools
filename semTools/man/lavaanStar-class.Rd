\name{lavaanStar-class}
\docType{class}
\alias{lavaanStar-class}
\alias{inspect,lavaanStar-method}
\alias{summary,lavaanStar-method}

\title{Class For Representing A (Fitted) Latent Variable Model with Additional Elements}
\description{This is the \code{lavaan} class that contains additional information about the fit values from the null model. Some functions are adjusted according to the change.}
\section{Objects from the Class}{
Objects can be created via the \code{\link{auxiliary}} function.
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
  }
}
% \section{Methods}{
  % \describe{
    % \item{summary}{\code{signature(object = "lavaan", type = "free")}: Returns 
      % the estimates of the parameters in the model as a named numeric vector.
      % If \code{type="free"}, only the free parameters are returned. If
      % \code{type="unco"}, both free and constrained parameters (simple
      % equality constraints only) are returned.
      % If \code{type="user"}, all parameters listed in the parameter table
      % are returned, including constrained and fixed parameters.}
    % \item{inspect}{\code{signature(object = "lavaan")}: Returns the 
      % implied moments of the model as a list with two elements (per group): 
      % \code{cov} for the implied covariance matrix, 
      % and \code{mean} for the implied mean 
      % vector. If only the covariance matrix was analyzed, the implied mean 
      % vector will be zero.}
  % }
% }
\references{
see \code{linkS4class{lavaan}}}
\seealso{
\code{\link{auxiliary}}
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
			  
dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
			  
fit <- cfa(HS.model, data=dat) #, group="sex", meanstructure=TRUE)
fitaux <- auxiliary(fit, aux="z", data=dat)
}
