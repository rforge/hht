\name{dcb_extractimf}
\alias{dcb_extractimf}
\title{Intrinsic Mode Function}
\description{
This function extracts intrinsic mode function from a signal. Most of this code is from the ``extractimf'' function in the EMD package.
}
\usage{
dcb_extractimf(residue, tt=NULL, tol=sd(residue)*0.1^2, max.sift=20, 
    stoprule="type1", boundary="periodic", sm="none", spar=NA, 
    weight=20, check=FALSE)
}
\arguments{
\item{residue}{observation or signal observed at time \code{tt}}
\item{tt}{observation index or time index}
\item{tol}{tolerance for stopping rule of sifting}
\item{max.sift}{the maximum number of sifting}
\item{stoprule}{stopping rule of sifting: ``type1'' (envelope mean <=tol), ``type2'' (standard deviation between two siftings <= tol), ``type3'' (S stoppage criterion)}
\item{boundary}{specifies boundary condition from ``none'', ``wave'', ``symmetric'', ``periodic'' or ``evenodd''.}
\item{sm}{specifies whether envelope is constructed by smoothing spline.} 
\item{spar}{specifies user-supplied smoothing parameter of spline.}
\item{weight}{the smoothness of spline is determined by \code{weight} times smoothing parameter of GCV.}
\item{check}{specifies whether the sifting process is displayed. If \code{check=TRUE}, click the plotting area to start the next step.}
}
\details{
This function extracts intrinsic mode functions from a signal.
}
\value{ 
\item{imf}{imf}
\item{residue}{residue signal after extracting the finest imf from \code{residue}}
\item{niter}{the number of iteration to obtain the \code{imf}} 
}
\references{
Huang, N. E., Shen, Z., Long, S. R., Wu, M. L. Shih, H. H.,
Zheng, Q., Yen, N. C., Tung, C. C. and Liu, H. H. (1998) The empirical
mode decomposition and Hilbert spectrum for nonlinear and
nonstationary time series analysis. \emph{Proceedings of the Royal
Society London A}, \bold{454}, 903--995.

Huang, N. E. and Wu Z. A. (2008) A review on Hilbert-Huang Transform: Method and its applications to geophysical studies.
\emph{Reviews of Geophysics}, \bold{46}, RG2006.
}

\author{Donghoh Kim \email{donghoh.kim@gmail.com}, Hee-Seok Oh, Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{dcb_emd}}
}
\examples{
### Generating a signal
ndata <- 3000
X11(); par(mfrow=c(1,1), mar=c(1,1,1,1))
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) + sin(6 * pi * tt2)  + 0.5 * tt2
plot(tt2, xt2, xlab="", ylab="", type="l", axes=FALSE); box()

### Extracting the first IMF by sifting process
tryimf <- dcb_extractimf(xt2, tt2, check=FALSE)
}
\keyword{nonparametric}
