\name{precision_tester}
\alias{precision_tester}
\title{
Test numerically determined instantaneous frequency against exact instantaneous frequency
}
\description{
This function compares the performance of \code{\link{instantaneous_frequency}} against signals of known instantaneous frequency.
The known signal is of the form:
\deqn{ x(t) = a\sin\{\omega_\{1\} + \varphi_\{1\}} + b\sin\{\omega_\{2\} + \varphi_\{2\}} + c}{a * sin(omega_1 + phi_1) + b * sin(omega_2 + phi_2) + c}
One can create quite complicated signals by choosing the various amplitude, frequency, and phase constants.

}
\usage{
precision_tester(tt = seq(0, 10, by = 0.01), method = "arctan", lag = 1, a = 1, b = 1, c = 1, omega_1 = 2 * pi, omega_2 = 4 * pi, phi_1 = 0, phi_2 = pi/6, plot_signal = TRUE, plot_instfreq = TRUE, plot_error = TRUE, ...)
}
\arguments{
  \item{tt}{
Sample times.
}
  \item{method}{
How the numeric instantaneous frequency is calculated, see \code{\link{instantaneous_frequency}}
}
  \item{lag}{
Differentiation lag, see the \code{diff} function in the \code{base} package.
}
  \item{a}{
Amplitude coefficient for the first sinusoid.
}
  \item{b}{
Amplitude coefficient for the second sinusoid.
}
  \item{c}{
DC shift
}
  \item{omega_1}{
Frequency of the first sinusoid.
}
  \item{omega_2}{
Frequency of the second sinusoid.
}
  \item{phi_1}{
Phase shift of the first sinusoid.
}
  \item{phi_2}{
Phase shift of the second sinusoid.
}
  \item{plot_signal}{
Whether to open a new window and show the time series.
}
  \item{plot_instfreq}{
Whether to open a new plot window and show the instantaneous frequencies, comparing the numerical and analytical result.
}
  \item{plot_error}{
Whether to open a new plot window and show the difference between the numerical and analytical result.
}
  \item{\dots}{
Plotting parameters.
}
}
\details{
For a description of how the exact analytical frequency is derived, see \link{http://www.unc.edu/\%7Ehaksaeng/hht/analytic_instantaneous_freq.pdf}
}
\value{
\item{instfreq$sig}{The time series}
\item{instfreq$analytic}{The exact instantaneous frequency}
\item{instfreq$numeric}{The numerically-derived instantaneous frequency from \code{\link{instantaneous_frequency}}}
}
\author{
Daniel C. Bowman \email{daniel.bowman@unc.edu}
}

\examples{

#Simple signal
tt = seq(0, 10, by = 0.01)
a = 1
b = 0
c = 0
omega_1 = 30 * pi
omega_2 = 0
phi_1 = 0
phi_2 = 0
precision_tester(tt, method = "arctan", lag = 1, a, b, c, omega_1, omega_2, phi_1, phi_2)

#That was nice - what happens if we use the "chain" method...?

precision_tester(tt, method = "chain", lag = 1, a, b, c, omega_1, omega_2, phi_1, phi_2)

#Big problems!  Let's increase the sample rate

tt = seq(0, 10, by = 0.0005)
precision_tester(tt, method = "chain", lag = 1, a, b, c, omega_1, omega_2, phi_1, phi_2)
#That's better

#Frequency modulations caused by signal that is not symmetric about 0

tt = seq(0, 10, by = 0.01)
a = 1
b = 0
c = 0.25
omega_1 = 2 * pi
omega_2 = 0
phi_1 = 0
phi_2 = 0

precision_tester(tt, method = "arctan", lag = 1, a, b, c, omega_1, omega_2, phi_1, phi_2)

#Non-uniform sample rate
set.seed(628)
tt = sort(runif(500, 0, 10))
a = 1
b = 0
c = 0
omega_1 = 2 * pi
omega_2 = 0
phi_1 = 0
phi_2 = 0

precision_tester(tt, method = "arctan", lag = 1, a, b, c, omega_1, omega_2, phi_1, phi_2)
}

\seealso{\code{\link{instantaneous_frequency}}}

\keyword{nonparametric}
