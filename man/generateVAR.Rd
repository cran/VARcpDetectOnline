% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_function.R
\name{generateVAR}
\alias{generateVAR}
\title{Generate VAR Data}
\usage{
generateVAR(n, As, Sig, h, isOldProvided = FALSE, oldxs = NULL)
}
\arguments{
\item{n}{Integer. The length of the VAR data to be generated.}

\item{As}{List. A list containing the transition matrices for the VAR process.}

\item{Sig}{Matrix. The covariance matrix of errors.}

\item{h}{Integer. The order of the VAR process.}

\item{isOldProvided}{Logical. If TRUE, the VAR data will be generated based on the
last observations from the previous segment of data. Defaults to FALSE.}

\item{oldxs}{Matrix. A \code{p} by \code{h} matrix containing the last observations from
the previous segment of data. Required if \code{isOldProvided = TRUE}.}
}
\value{
A data matrix of dimensions \code{p} by \code{n}.
}
\description{
This function generates Vector Auto-Regressive (VAR) data based on the provided parameters.
}
\examples{
# Example usage
As <- list(matrix(c(0.5, 0.2, 0.1, 0.4), 2, 2))
Sig <- diag(2)
data <- generateVAR(n = 100, As = As, Sig = Sig, h = 1)
}
