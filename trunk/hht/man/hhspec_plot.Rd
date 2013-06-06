\name{hhspec_plot}
\alias{hhspec_plot}
\title{Display Hilbert periodogram}
\description{
This function displays the Hilbert periodogram, with options to plot individual IMFs and also the Fourier periodogram for comparison.
}
\usage{
hhspec_plot(hspec, freq_span = NULL, scaling = "none", imf_list = NULL, 
show_total = TRUE, show_fourier = FALSE, show_imfs = FALSE, legend = TRUE, ...)
}
\arguments{
\item{hspec}{Data structure returned by \code{\link{hh_spectrum}}}
\item{freq_span}{Frequency range to plot, \code{NULL} plots all of them}
\item{scaling}{Amplitude scaling, can be \code{"log"} (log 10), \code{"sqrt"} (square root), defaults to \code{"none"}.}
\item{imf_list}{Which IMFs to plot, requires \code{show_imfs = TRUE}.}
\item{show_total}{Show the ensemble Hilbert spectrogram}
\item{show_fourier}{Show the Fourier periodogram}
\item{scale_fourier}{Scale Fourier and Hilbert spectra to each other for easier comparison}
\item{show_imfs}{Plot individual IMF spectra}
\item{legend}{Determines whether or not a legend is shown}
\item{This function supports some optional parameters as well:}
\itemize{
    \item{xlab}{X axis label}
    \item{ylab}{Y axis label}
    \item{legend_location}{where to put the legend}
    \item{total_col}{color of ensemble Hilbert periodogram}
    \item{total_lwd}{lwd of ensemble Hilbert periodogram}
    \item{total_lty}{lty of ensemble Hilbert periodogram}
    \item{imf_cols}{colors of IMF periodogram}
    \item{imf_lwd}{lwds of IMF periodogram}
    \item{imf_lty}{ltys of IMF periodogram}
    \item{fourier_col}{color of Fourier periodogram}
    \item{fourier_lwd}{lwd of Fourier periodogram}
    \item{fourier_lty}{lty of Fourier periodogram}
    \item{main}{Adds a title to the figure}}
}

\details{
This function plots the Hilbert periodogram of a signal, with options to show periodograms of individual IMFs.
You can also plot a simple Fourier periodogram for comparison.
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{hh_spectrum}}, \code{\link{hhgram_image}}}
\examples{

#Here we see how the EMD produces a dyadic filter bank for uniform random noise
#The frequency distributions of all but the first IMF display a Chi-Square distribution
#See Huang, N. E. & Wu, Z. 
#A review on Hilbert-Huang Transform: Method and its applications to geophysical studies
#Reviews of Geophysics, 2008, 46, RG2006

#The EMD of this signal may take a couple of minutes to run

set.seed(628)
sig = runif(10000)
tt = seq(length(sig)) * 0.01

\dontrun{emd_result = sig2imf(sig, tt)}

dfreq = 0.1
\dontrun{hspec = hh_spectrum(emd_result, dfreq)}

\dontrun{hhspec_plot(hspec, show_imfs = TRUE, 
imf_list = 1:10, show_total = TRUE, scaling = "sqrt", 
imf_lwd = rep(2, 10), total_lty = 3)}
}

\keyword{nonparametric}
