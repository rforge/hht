\name{combine_trials}
\alias{combine_trials}
\title{Gather trial files}
\description{
This function gathers trial files from multiple directories, renumbers them, and saves them to a single directory for processing using \code{\link{eemd_compile}}.
}
\usage{
combine_trials(in_dirs, out_dir, copy=TRUE)  
}
\arguments{
\item{in_dirs}{Directories containing trial file sets from one EEMD run.}
\item{out_dir}{Directory in which to save all trial files.}
\item{copy}{Copy files (copy=TRUE) or move them (copy=FALSE).}
}
\details{
Parallel processing is an efficient method for running EEMD.
However, this will result in several directories, each with trial files numbered from 1 to N.
These files cannot simply be copied together into the same directory, because then they would overwrite each other.
This function gathers all trial files in multiple directories, renumbers them, and saves them in a different directory.
}
\value{ 
The trial files are saved in the directory specified by \code{out_dir}.
}

\author{Daniel Bowman \email{daniel.bowman@unc.edu}}

\seealso{
\code{\link{eemd}}, \code{\link{eemd_compile}}
}
\examples{
#Suppose you have run 3 different EEMD sets of 100 trials each and saved the results in eemd1, eemd2, eemd3, respectively:
in_dirs=c("/home/user/eemd1", "/home/user/eemd2/", "/home/user/eemd3")
out_dir="/home/user/all_trials"
\dontrun{combine_trials(in_dirs, out_dir)}
#Now all your trials should be located in /home/user/all_trials, numbered 1 through 300
}
\keyword{nonparametric}
