\name{DecomposeResidue}
\alias{DecomposeResidue}
\title{
Extract IMFs from residue
}
\description{
Ensure that residues returned by \code{\link{Sig2IMF}} have at most one extremum.
This is done by differentiating the residue produced by the EMD method, running the EMD again on the differentiated remainder, then reintegrating.
}

\usage{
DecomposeResidue(emd.result)
}
\arguments{
  \item{emd.result}{
  Output from \code{\link{Sig2IMF}}}
}
\value{
\item{emd.result}{The result of the EMD, now with a residue that has at most one extremum.}
}
\author{
Daniel Bowman \email{daniel.bowman@unc.edu}
}
\note{
The numerical differentiation and integration will introduce error, meaning that the sum of the new IMF set and the residual will not exactly match the original signal.
This function is only used internally and will not be called by the user under normal circumstances.}
