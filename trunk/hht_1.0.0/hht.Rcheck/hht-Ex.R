pkgname <- "hht"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('hht')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("combine_trials")
### * combine_trials

flush(stderr()); flush(stdout())

### Name: combine_trials
### Title: Gather trial files
### Aliases: combine_trials
### Keywords: nonparametric

### ** Examples

#Suppose you have run 3 different EEMD sets of 100 trials each and saved the results in eemd1, eemd2, eemd3, respectively:
in_dirs=c("/home/user/eemd1", "/home/user/eemd2/", "/home/user/eemd3")
out_dir="/home/user/all_trials"
## Not run: combine_trials(in_dirs, out_dir)
#Now all your trials should be located in /home/user/all_trials, numbered 1 through 300



cleanEx()
nameEx("dcb_emd")
### * dcb_emd

flush(stderr()); flush(stdout())

### Name: dcb_emd
### Title: Empirical Mode Decomposition
### Aliases: dcb_emd
### Keywords: nonparametric

### ** Examples

### Empirical Mode Decomposition
ndata <- 3000
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) + sin(6 * pi * tt2)  + 0.5 * tt2

par(mfrow=c(3,1), mar=c(2,1,2,1))
try <- dcb_emd(xt2, tt2, boundary="wave")

### Plotting the IMFs
par(mfrow=c(3,1), mar=c(2,1,2,1))
X11(); par(mfrow=c(try$nimf+1, 1), mar=c(2,1,2,1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf) {
    plot(tt2, try$imf[,i], type="l", xlab="", ylab="", ylim=rangeimf,
    main=paste(i, "-th IMF", sep="")); abline(h=0)
}
plot(tt2, try$residue, xlab="", ylab="", main="residue", type="l", 
axes=FALSE); box()



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("dcb_extractimf")
### * dcb_extractimf

flush(stderr()); flush(stdout())

### Name: dcb_extractimf
### Title: Intrinsic Mode Function
### Aliases: dcb_extractimf
### Keywords: nonparametric

### ** Examples

### Generating a signal
ndata <- 3000
X11(); par(mfrow=c(1,1), mar=c(1,1,1,1))
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) + sin(6 * pi * tt2)  + 0.5 * tt2
plot(tt2, xt2, xlab="", ylab="", type="l", axes=FALSE); box()

### Extracting the first IMF by sifting process
tryimf <- dcb_extractimf(xt2, tt2, check=FALSE)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("eemd")
### * eemd

flush(stderr()); flush(stdout())

### Name: eemd
### Title: Ensemble Empirical Mode Decomposition
### Aliases: eemd
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_config$S=5

trials=10
nimf=10
noise_amp=6.4e-07
trials_dir="test"

set.seed(628)
#Run EEMD (this may take some time)
## Not run: eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)

#Compile the results
## Not run: eemd_result=eemd_compile(trials_dir, trials, nimf)

#Plot the IMFs
time_span=c(5, 10)
imf_list=1:3
os=TRUE
res=TRUE
## Not run: plot_imfs(eemd_result, time_span, imf_list, os, res)



cleanEx()
nameEx("eemd_compile")
### * eemd_compile

flush(stderr()); flush(stdout())

### Name: eemd_compile
### Title: Process EEMD results
### Aliases: eemd_compile
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_config$S=5

trials=10
nimf=10
noise_amp=6.4e-07
trials_dir="test"
set.seed(628)
#Run EEMD (this may take some time)
## Not run: eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)

#Compile the results
## Not run: eemd_result=eemd_compile(trials_dir, trials, nimf)

#Plot the IMFs
time_span=c(5, 10)
imf_list=1:3
os=TRUE
res=TRUE
## Not run: plot_imfs(eemd_result, time_span, imf_list, os, res)



cleanEx()
nameEx("eemd_resift")
### * eemd_resift

flush(stderr()); flush(stdout())

### Name: eemd_resift
### Title: Resift averaged IMFs from EEMD
### Aliases: eemd_resift
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20

trials=10
nimf=10
noise_amp=6.4e-07
trials_dir="test"

set.seed(628)

#Run EEMD (this may take some time)
## Not run: eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)

#Compile the results
## Not run: eemd_result=eemd_compile(trials_dir, trials, nimf)


resift_rule="max_var"
## Not run: resift_result=eemd_resift(eemd_result, emd_config, resift_rule)

#Plot the IMFs
time_span=c(5, 10)
imf_list=1:3
os=TRUE
res=TRUE
## Not run: plot_imfs(resift_result, time_span, imf_list, os, res)



cleanEx()
nameEx("evolutive_fft")
### * evolutive_fft

flush(stderr()); flush(stdout())

### Name: evolutive_fft
### Title: Generate Fourier spectrogram
### Aliases: evolutive_fft
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

ft=list()
ft$nfft=4096
ft$ns=30
ft$nov=29

freq_span=c(0, 25)
ev = evolutive_fft(sig, dt, ft, freq_span)

#Plot raw spectrogram
image_z=t(ev$DSPEC[1:(ev$numfreqs/2),])
f_ind=(ev$freqs>=freq_span[1] & ev$freqs <= freq_span[2])
image_z=t(ev$DSPEC[1:(ev$numfreqs/2),])
image_z=image_z[,f_ind]
image(image_z)



cleanEx()
nameEx("ftspec_image")
### * ftspec_image

flush(stderr()); flush(stdout())

### Name: ftspec_image
### Title: Display Fourier spectrogram
### Aliases: ftspec_image
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

ft=list()
ft$nfft=4096
ft$ns=30
ft$nov=29

time_span=c(5, 10)
freq_span=c(0, 25)
amp_span=c(1e-5, 0.0003)
ftspec_image(sig, dt, ft, time_span, freq_span, amp_span)



cleanEx()
nameEx("hh_render")
### * hh_render

flush(stderr()); flush(stdout())

### Name: hh_render
### Title: Render Hilbert spectrogram
### Aliases: hh_render
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_config$S=5

trials=10
nimf=10
noise_amp=6.4e-07
trials_dir="test"

set.seed(628)
#Run EEMD (this may take some time)
## Not run: eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)

#Compile the results
## Not run: eemd_result=eemd_compile(trials_dir, trials, nimf)

#Calculate spectrogram
max_freq=25
freq_step=0.01
## Not run: hspec=hh_render(eemd_result, max_freq, freq_step)

#Plot spectrogram 
time_span=c(5, 10)
freq_span=c(0, 25)
amp_span=c(1e-6, 2.5e-5)
## Not run: hhspec_image(hspec, time_span, freq_span, amp_span)



cleanEx()
nameEx("hhspec_image")
### * hhspec_image

flush(stderr()); flush(stdout())

### Name: hhspec_image
### Title: Display Hilbert Huang spectrogram
### Aliases: hhspec_image
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_config$S=5

trials=10
nimf=10
noise_amp=6.4e-07
trials_dir="test"

set.seed(628)
#Run EEMD (this may take some time)
## Not run: eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)

#Compile the results
## Not run: eemd_result=eemd_compile(trials_dir, trials, nimf)

#Calculate spectrogram
max_freq=25
freq_step=0.01
## Not run: hspec=hh_render(eemd_result, max_freq, freq_step)

#Plot spectrogram 
time_span=c(5, 10)
freq_span=c(0, 25)
amp_span=c(1e-6, 2.5e-5)
## Not run: hhspec_image(hspec, time_span, freq_span, amp_span)



cleanEx()
nameEx("hhtransform")
### * hhtransform

flush(stderr()); flush(stdout())

### Name: hhtransform
### Title: Hilbert transform wrapper
### Aliases: hhtransform
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_config$S=5

#Run EMD (this may take some time)
emd_result=sig2imf(sig, dt, emd_config)

#Get instantaneous amplitude and frequency
emd_result=hhtransform(emd_result)

#Render spectrogram
max_freq=25
freq_step=0.05
hspec=hh_render(emd_result, max_freq, freq_step)

#Show result
time_span=c(5, 10)
freq_span=c(0, 25)
amp_span=c(0.000001, 0.00001)
hhspec_image(hspec, time_span, freq_span, amp_span)



cleanEx()
nameEx("plot_imfs")
### * plot_imfs

flush(stderr()); flush(stdout())

### Name: plot_imfs
### Title: Display IMFs
### Aliases: plot_imfs
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=5
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20

#Run EMD
emd_result=sig2imf(sig, dt, emd_config)

#Plot the first 4 IMFs of the EEMD of a signal.
time_span=c(5, 10)
imf_list=1:4
original_signal=TRUE
residue=TRUE

plot_imfs(emd_result, time_span, imf_list, original_signal, residue)

#Check how much contribution IMFs 2 and 3 make to the complete signal.
imf_list=c(2, 3)
fit_line=TRUE
plot_imfs(emd_result, time_span, imf_list, original_signal, residue, fit_line)



cleanEx()
nameEx("sig2imf")
### * sig2imf

flush(stderr()); flush(stdout())

### Name: sig2imf
### Title: Empirical Mode Decomposition wrapper
### Aliases: sig2imf
### Keywords: nonparametric

### ** Examples

data(port_foster_event)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=5
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20

#Run EMD
emd_result=sig2imf(sig, dt, emd_config)

#Display IMFs

time_span=c(5, 10)
imf_list=1:3
original_signal=TRUE
residue=TRUE

plot_imfs(emd_result, time_span, imf_list, original_signal, residue)

#Get Hilbert transform
emd_result=hhtransform(emd_result)

#Plot spectrogram
max_freq=25
freq_step=0.05
hspec=hh_render(emd_result, max_freq, freq_step)

freq_span=c(0, 25)
amp_span=c(0.000001, 0.00001)
hhspec_image(hspec, time_span, freq_span, amp_span)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
