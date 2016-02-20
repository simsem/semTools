\name{indProd}
\alias{indProd}
\alias{orthogonalize}
\title{
	Make products of indicators using no centering, mean centering, double-mean centering, or residual centering
}
\description{
The \code{indProd} function will make products of indicators using no centering, mean centering, double-mean centering, or residual centering. The \code{orthogonalize} function is the shortcut of the \code{indProd} function to make the residual-centered indicators products.
}
\usage{
indProd(data, var1, var2, var3=NULL, match = TRUE, meanC = TRUE, 
	residualC = FALSE, doubleMC = TRUE, namesProd = NULL)
orthogonalize(data, var1, var2, var3=NULL, match=TRUE, namesProd=NULL)
}
\arguments{
  \item{data}{
	The desired data to be transformed.
}
  \item{var1}{
	Names or indices of the variables loaded on the first factor
}
  \item{var2}{
	Names or indices of the variables loaded on the second factor
}
  \item{var3}{
	Names or indices of the variables loaded on the third factor (for three-way interaction)
}
  \item{match}{
	Specify \code{TRUE} to use match-paired approach (Marsh, Wen, & Hau, 2004). If \code{FALSE}, the resulting products are all possible products.
}
  \item{meanC}{
	Specify \code{TRUE} for mean centering the main effect indicator before making the products
}
  \item{residualC}{
	Specify \code{TRUE} for residual centering the products by the main effect indicators (Little, Bovaird, & Widaman, 2006).
}
  \item{doubleMC}{
	Specify \code{TRUE} for centering the resulting products (Lin et. al., 2010)
}
  \item{namesProd}{
	The names of resulting products
}
}
\value{
	The original data attached with the products.
}
\references{
	Marsh, H. W., Wen, Z. & Hau, K. T. (2004). Structural equation models of latent interactions: Evaluation of alternative estimation strategies and indicator construction. \emph{Psychological Methods, 9,} 275-300.
	
	Lin, G. C., Wen, Z., Marsh, H. W., & Lin, H. S. (2010). Structural equation models of latent interactions: Clarification of orthogonalizing and double-mean-centering strategies. \emph{Structural Equation Modeling, 17}, 374-391.
	
	Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of orthogonalizing powered and product terms: Implications for modeling interactions among latent variables. \emph{Structural Equation Modeling, 13}, 497-519.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
	Alexander Schoemann (East Carolina University; \email{schoemanna@ecu.edu})
}
\seealso{
	\itemize{
		\item \code{\link{probe2WayMC}} For probing the two-way latent interaction when the results are obtained from mean-centering, or double-mean centering.
		\item \code{\link{probe3WayMC}} For probing the three-way latent interaction when the results are obtained from mean-centering, or double-mean centering.
		\item \code{\link{probe2WayRC}} For probing the two-way latent interaction when the results are obtained from residual-centering approach.
		\item \code{\link{probe3WayRC}} For probing the two-way latent interaction when the results are obtained from residual-centering approach.
		\item \code{\link{plotProbe}} Plot the simple intercepts and slopes of the latent interaction.
	}
}
\examples{
# Mean centering / two-way interaction / match-paired
dat <- indProd(attitude[,-1], var1=1:3, var2=4:6)

# Residual centering / two-way interaction / match-paired
dat2 <- indProd(attitude[,-1], var1=1:3, var2=4:6, match=FALSE, meanC=FALSE, 
	residualC=TRUE, doubleMC=FALSE)

# Double-mean centering / two-way interaction / match-paired
dat3 <- indProd(attitude[,-1], var1=1:3, var2=4:6, match=FALSE, meanC=TRUE, 
	residualC=FALSE, doubleMC=TRUE)

# Mean centering / three-way interaction / match-paired
dat4 <- indProd(attitude[,-1], var1=1:2, var2=3:4, var3=5:6)

# Residual centering / three-way interaction / match-paired
dat5 <- indProd(attitude[,-1], var1=1:2, var2=3:4, var3=5:6, match=FALSE, meanC=FALSE, 
	residualC=TRUE, doubleMC=FALSE)

# Double-mean centering / three-way interaction / match-paired
dat6 <- indProd(attitude[,-1], var1=1:2, var2=3:4, var3=5:6, match=FALSE, meanC=TRUE, 
	residualC=TRUE, doubleMC=TRUE)
}
