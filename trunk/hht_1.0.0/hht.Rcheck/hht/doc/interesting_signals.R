### R code from vignette source 'interesting_signals.Rnw'

###################################################
### code chunk number 1: interesting_signals.Rnw:48-64 (eval = FALSE)
###################################################
## library(hht)
## dt=0.01
## tt=seq_len(10000)*dt
## sig=0.5*sin(pi*tt/10)+sin(2*pi*tt)+0.5*sin(10*pi*tt)+0.5*sin(40*pi*tt)
## emd_config=list()
## emd_config$max_sift=200
## emd_config$max_imf=100
## emd_config$tol=5
## emd_config$stop_rule="type5"
## emd_config$boundary="wave"
## emd_config$sm="none"
## emd_config$spar=NA
## emd_config$weight=20
## 
## emd_result=sig2imf(sig, dt, emd_config)
## emd_result=hhtransform(emd_result)


###################################################
### code chunk number 2: interesting_signals.Rnw:68-73 (eval = FALSE)
###################################################
## time_span=c(0,40)
## imf_list=1:4
## os=TRUE
## res=FALSE
## plot_imfs(emd_result, time_span, imf_list, os, res)


###################################################
### code chunk number 3: interesting_signals.Rnw:78-81 (eval = FALSE)
###################################################
## max_freq=25
## freq_step=0.1
## hspec=hh_render(emd_result, max_freq, freq_step)


###################################################
### code chunk number 4: interesting_signals.Rnw:84-88 (eval = FALSE)
###################################################
## time_span=c(0, -1)
## freq_span=c(0, 25)
## amp_span=c(0, 1)
## hhspec_image(hspec, time_span, freq_span, amp_span, cex=1)


###################################################
### code chunk number 5: interesting_signals.Rnw:93-99 (eval = FALSE)
###################################################
## ft=list()
## ft$nfft=2048
## ft$ns=1000
## ft$nov=998
## amp_span=c(100,200)
## ftspec_image(sig, dt, ft, time_span, freq_span, amp_span, cex=1)


###################################################
### code chunk number 6: interesting_signals.Rnw:154-174 (eval = FALSE)
###################################################
## library(hht)
## dt=0.05
## tt=seq_len(10000)*dt
## sig=sin(2*pi*tt)+0.25*sin(4*pi*tt)
## emd_config=list()
## emd_config$max_sift=200
## emd_config$max_imf=100
## emd_config$tol=5
## emd_config$stop_rule="type5"
## emd_config$boundary="wave"
## emd_config$sm="none"
## emd_config$spar=NA
## emd_config$weight=20
## emd_result=sig2imf(sig, dt, emd_config)
## time_span=c(20, 30)
## imf_list=c(1)
## os=TRUE
## res=FALSE
## fit_line=TRUE
## plot_imfs(emd_result, time_span, imf_list, os, res, fit_line)


###################################################
### code chunk number 7: interesting_signals.Rnw:179-184 (eval = FALSE)
###################################################
## emd_result=hhtransform(emd_result)
## max_freq=4
## freq_step=0.1
## imf_list=c(1)
## hspec=hh_render(emd_result, max_freq, freq_step, imf_list)


###################################################
### code chunk number 8: interesting_signals.Rnw:187-191 (eval = FALSE)
###################################################
## time_span=c(20, 30)
## freq_span=c(0, 3)
## amp_span=c(0, -1)
## hhspec_image(hspec, time_span, freq_span, amp_span, cex=1)


###################################################
### code chunk number 9: interesting_signals.Rnw:195-203 (eval = FALSE)
###################################################
## ft=list()
## ft$nfft=4096
## ft$ns=2000
## ft$nov=1998
## time_span=c(0, -1)
## freq_span=c(0, 3)
## amp_span=c(90, 200)
## ftspec_image(sig, dt, ft, time_span, freq_span, amp_span, cex=1)


###################################################
### code chunk number 10: nonlinstaimfs
###################################################
library(hht)
dt=0.01
tt=seq_len(10000)*dt
sig=sin(2*pi*tt+0.5*sin(pi*tt))*exp(-(tt-50)^2/100)+sin(pi*tt/10)
emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=5
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_result=sig2imf(sig, dt, emd_config)

time_span=c(0, -1)
imf_list=c(1:2)
os=TRUE
res=TRUE
plot_imfs(emd_result, time_span, imf_list, os, res)


###################################################
### code chunk number 11: interesting_signals.Rnw:286-290
###################################################
emd_result=hhtransform(emd_result)
max_freq=5
freq_step=0.05
hspec=hh_render(emd_result, max_freq, freq_step)


###################################################
### code chunk number 12: nonlinstaht
###################################################
time_span=c(0, 100)
freq_span=c(0, 2.5)
amp_span=c(0, 1)
hhspec_image(hspec, time_span, freq_span, amp_span, cex=1)


###################################################
### code chunk number 13: nonlinstaft
###################################################
ft=list()
ft$nfft=4096
ft$ns=2000
ft$nov=1998
amp_span=c(100,200)
ftspec_image(sig, dt, ft, time_span, freq_span, amp_span, cex=1)


###################################################
### code chunk number 14: nonlinstaimfs
###################################################
library(hht)
dt=0.01
tt=seq_len(10000)*dt
sig=sin(2*pi*tt+0.5*sin(pi*tt))*exp(-(tt-50)^2/100)+sin(pi*tt/10)
emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=5
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_result=sig2imf(sig, dt, emd_config)

time_span=c(0, -1)
imf_list=c(1:2)
os=TRUE
res=TRUE
plot_imfs(emd_result, time_span, imf_list, os, res)


###################################################
### code chunk number 15: nonlinstaht
###################################################
time_span=c(0, 100)
freq_span=c(0, 2.5)
amp_span=c(0, 1)
hhspec_image(hspec, time_span, freq_span, amp_span, cex=1)


###################################################
### code chunk number 16: nonlinstaft
###################################################
ft=list()
ft$nfft=4096
ft$ns=2000
ft$nov=1998
amp_span=c(100,200)
ftspec_image(sig, dt, ft, time_span, freq_span, amp_span, cex=1)


