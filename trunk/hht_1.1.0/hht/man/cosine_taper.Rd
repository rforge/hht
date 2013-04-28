\name{cosine_taper}
\alias{cosine_taper}
\title{Taper signal with cosine function}
\description{
This function applies a cosine taper to a time series.
It is a copy of the \code{rsspec.taper} function in the RSEIS package.
}
\usage{
cosine_taper(sig, taper = 0.05) 
}
\arguments{
\item{sig}{The signal to process}
\item{taper}{Fraction of signal to taper (greater than 0, less than 0.5)}
}

\references{Jonathan M. Lees (2012). RSEIS: Seismic Time Series Analysis Tools. R
  package version 3.0-6. http://CRAN.R-project.org/package=RSEIS}

\author{Jonathan M. Lees}

\examples{
dt = 0.01
tt = 1:1000 * dt
sig = sin(4 * pi * tt)
tapsig = cosine_taper(sig, taper = 0.25)
plot(tt, sig, type = "l")
points(tt, tapsig, col = "red")
}
\keyword{nonparametric}
