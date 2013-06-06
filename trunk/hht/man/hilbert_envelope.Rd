\name{hilbert_envelope}
\alias{hilbert_envelope}
\title{
Instantaneous amplitude
}
\description{
Generates the instantaneous amplitude of an analytic signal given by \code{\link{hilbert_transform}}
}
\usage{
hilbert_envelope(asig)
}
\arguments{
  \item{asig}{
The analytic signal returned by \code{\link{hilbert_transform}}
}
}
\value{
\item{envelope}{Instantaneous amplitude}
}
\author{
Daniel C. Bowman \email{daniel.bowman@unc.edu}
}

\seealso{
\code{\link{hilbert_transform}}, \code{\link{instantaneous_frequency}}
}
\examples{
tt = seq(1000) * 0.01
sig = sin(4 * pi * tt) + sin(3.4 * pi * tt)
asig = hilbert_transform(sig)
env = hilbert_envelope(asig)
plot(tt, sig, type = "l")
lines(tt, env, col = "red")
lines(tt, -env, col = "red")

}
\keyword{nonparametric}
