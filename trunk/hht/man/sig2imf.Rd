\name{sig2imf}
\alias{sig2imf}
\title{Empirical Mode Decomposition wrapper}
\description{
This function wraps the \code{emd} function in the \code{EMD} package.
\code{sig2imf} is used in \code{\link{eemd}} and others.
}
\usage{
sig2imf(sig, tt, spectral_method = "arctan", diff_lag = 1, stop_rule = "type5", 
    tol = 5, boundary = "wave", sm = "none", smlevels = c(1), spar = NULL, 
    max_sift = 200, max_imf = 100, interm = NULL)
}

\arguments{
\item{sig}{a time series to be decomposed (vector)}
\item{tt}{A vector of sample times for \code{sig}}
\item{spectral_method}{defines how to calculate instantaneous frequency - whether to use the arctangent of the analytic signal with numeric differentiation (``arctan'')
    or the result of the chain rule applied to the arctangent, then numerically differentiated ("chain");  see \code{\link{instantaneous_frequency}}.}
\item{diff_lag}{specifies if you want to do naive differentiation (DIFF_LAG = 1), central difference method (DIFF_LAG = 2) or higher difference methods (DIFF_LAG > 2) 
     to determine instantaneous frequency; see \code{\link{instantaneous_frequency}}.}
\item{stop_rule}{As quoted from the EMD package documentation:  ''The stop rule of sifting. 
     The type1 stop rule indicates that absolute values 
     of envelope mean must be less than the user-specified tolerance level in the sense
     that the local average of upper and lower envelope is zero. The stopping rules
     type2, type3, type4 and type5 are the stopping rules given by equation (5.5)
     of Huang et al. (1998), equation (11a), equation (11b) and S stoppage of Huang
     and Wu (2008), respectively.''}
\item{tol}{Determines what value is used to stop the sifting - this will depend on which stop rule you use.}
\item{boundary}{how the beginning and end of the signal are handled}
\item{sm}{Specifies how the signal envelope is constructed, see Kim et al, 2012.}
\item{smlevels}{Specifies what level of the IMF is obtained by smoothing other than interpolation, see EMD package documentation}
\item{spar}{User-defined smoothing parameter for spline, kernal, or local polynomial smoothing.}
\item{max_sift}{How many sifts are allowed - if this value is exceeded the IMF is returned as-is.}
\item{max_imf}{Maximum number of IMFs allowed.}
\item{interm}{Specifies vector of periods to be excluded from IMFs to cope with mode mixing.}
}
\details{
This function configures and performs empirical mode decomposition using the \code{emd} function in the \code{EMD} package.
}
\value{ 
\item{emd_result}{The instrinsic mode functions (IMFs), instantaneous frequencies, and instantaneous ampltitudes of \code{sig}.}
}
\references{
Kim, D., Kim, K. and Oh, H.-S. (2012) Extending the scope of empirical mode decomposition by smoothing.
\emph{EURASIP Journal on Advances in Signal Processing}, \bold{2012}, 168.

Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.

Huang, N. E. and Wu Z. A. (2008) A review on Hilbert-Huang Transform: Method and its applications to geophysical studies.
\emph{Reviews of Geophysics}, \bold{46}, RG2006.
}

\seealso{
\code{\link{eemd}}, \code{\link{plot_imfs}}
}
\examples{
data(port_foster_event)

#Run EMD
emd_result=sig2imf(sig, tt)

#Display IMFs

time_span <- c(5, 10)
imf_list <- 1:3
original_signal <- TRUE
residue <- TRUE

plot_imfs(emd_result, time_span, imf_list, original_signal, residue)

#Plot spectrogram
sdt <- tt[2] - tt[1]
dfreq <- 0.25
freq_span <- c(0, 25)
hgram <- hh_render(emd_result, sdt, dfreq, freq_span = freq_span, verbose = FALSE)

time_span <- c(4, 10)
freq_span <- c(0, 25)
amp_span <- c(0.000001, 0.00001)
hhgram_image(hgram, time_span = time_span, 
freq_span = freq_span, amp_span = amp_span,
pretty = TRUE)
}
\keyword{ts}
