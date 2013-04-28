### R code from vignette source 'eemd_examples.Rnw'

###################################################
### code chunk number 1: eemd_examples.Rnw:59-79 (eval = FALSE)
###################################################
## library(hht)
## dt=0.01
## tt=seq_len(10000)*dt
## sig=sin(2*pi*tt+0.5*sin(pi*tt))*exp(-(tt-50)^2/100)+sin(pi*tt/10)
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
## 
## time_span=c(0, -1)
## imf_list=c(1:2)
## os=TRUE
## res=TRUE
## plot_imfs(emd_result, time_span, imf_list, os, res)


###################################################
### code chunk number 2: eemd_examples.Rnw:83-92 (eval = FALSE)
###################################################
## set.seed(628)
## sig2=sin(2*pi*tt+0.5*sin(pi*tt))*exp(-(tt-50)^2/100)+
##     sin(pi*tt/10)+rnorm(length(sig), 0, 0.01)
## emd_result=sig2imf(sig2, dt, emd_config)
## time_span=c(0, -1)
## imf_list=c(1:emd_result$nimf)
## os=TRUE
## res=FALSE
## plot_imfs(emd_result, time_span, imf_list, os, res)


###################################################
### code chunk number 3: eemd_examples.Rnw:98-105 (eval = FALSE)
###################################################
## trials=100
## nimf=11
## noise_amp=0.01
## trials_dir="example_1"
##  
## eemd(sig2, dt, trials, nimf, noise_amp, emd_config, trials_dir)
## eemd_result=eemd_compile(trials_dir, trials, nimf)


###################################################
### code chunk number 4: eemd_examples.Rnw:110-128 (eval = FALSE)
###################################################
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
## resift_rule="max_var"
## resift_result=eemd_resift(eemd_result, emd_config, resift_rule)
## 
## time_span=c(0, 100)
## imf_list=c(3:7)
## os=TRUE
## res=FALSE
## plot_imfs(resift_result, time_span, imf_list, os, res)


###################################################
### code chunk number 5: eemd_examples.Rnw:196-216 (eval = FALSE)
###################################################
## library(hht)
## set.seed(628)
## dt=0.01
## tt=seq_len(2000)*dt
## sig=sin(2*pi*tt) + 2*sin(12*pi*tt)+rnorm(length(tt), 0, 0.4)
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
## time_span=c(0, 20)
## imf_list=1:5
## os=TRUE
## res=FALSE
## plot_imfs(emd_result, time_span, imf_list, os, res)


###################################################
### code chunk number 6: eemd_examples.Rnw:221-230 (eval = FALSE)
###################################################
## emd_result=hhtransform(emd_result)
## max_freq=10
## freq_step=0.1
## hspec=hh_render(emd_result, max_freq, freq_step)
## 
## time_span=c(0, 20)
## freq_span=c(0, 10)
## amp_span=c(0.75, 3)
## hhspec_image(hspec, time_span, freq_span, amp_span)


###################################################
### code chunk number 7: eemd_examples.Rnw:235-242 (eval = FALSE)
###################################################
## trials=100
## nimf=9
## noise_amp=0.4
## trials_dir="example_2"
## 
## eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)
## eemd_result=eemd_compile(trials_dir, trials, nimf)


###################################################
### code chunk number 8: eemd_examples.Rnw:247-272 (eval = FALSE)
###################################################
## 
## trials=100
## nimf=9
## trials_dir="example_2"
## eemd_result=eemd_compile(trials_dir, trials, nimf)
## 
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
## resift_rule="max_var"
## resift_result=eemd_resift(eemd_result, emd_config, resift_rule)
## 
## time_span=c(0, 20)
## imf_list=1:5
## os=TRUE
## res=FALSE
## 
## plot_imfs(resift_result, time_span, imf_list, os, res)


###################################################
### code chunk number 9: eemd_examples.Rnw:277-291 (eval = FALSE)
###################################################
## trials_dir="example_2"
## trials=100
## nimf=9
## eemd_result=eemd_compile(trials_dir, trials, nimf)
## 
## max_freq=10
## freq_step=0.1
## hspec=hh_render(eemd_result, max_freq, freq_step)
## 
## time_span=c(0, 20)
## freq_span=c(0, 10)
## amp_span=c(0.75, 3)
## cluster_span=c(8, 100)
## hhspec_image(hspec, time_span, freq_span, amp_span, cluster_span)


###################################################
### code chunk number 10: eemd_examples.Rnw:359-378 (eval = FALSE)
###################################################
## library(hht)
## data(port_foster_event)
## set.seed(628)
## emd_config=list()
## emd_config$max_sift=200
## emd_config$max_imf=100
## emd_config$tol=5
## emd_config$stop_rule="type5"
## emd_config$boundary="wave"
## emd_config$sm="none"
## emd_config$spar=NA
## emd_config$weight=20
## trials=100
## nimf=8
## noise_amp=6.4e-07
## trials_dir="example_3"
## 
## eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)
## eemd_result=eemd_compile(trials_dir, trials, nimf)


###################################################
### code chunk number 11: eemd_examples.Rnw:383-400 (eval = FALSE)
###################################################
## emd_config=list()
## emd_config$max_sift=200
## emd_config$max_imf=100
## emd_config$tol=5
## emd_config$stop_rule="type5"
## emd_config$boundary="wave"
## emd_config$sm="none"
## emd_config$spar=NA
## emd_config$weight=20
## resift_rule="max_var"
## resift_result=eemd_resift(eemd_result, emd_config, resift_rule)
## time_span=c(5, 10)
## imf_list=1:3
## os=TRUE
## res=FALSE
## fit_line=TRUE
## plot_imfs(resift_result, time_span, imf_list, os, res, fit_line)


###################################################
### code chunk number 12: eemd_examples.Rnw:404-419 (eval = FALSE)
###################################################
## trials_dir="example_3"
## trials=100
## nimf=11
## 
## eemd_result=eemd_compile(trials_dir, trials, nimf)
## 
## max_freq=30
## freq_step=0.05
## hspec=hh_render(eemd_result, max_freq, freq_step)
## 
## time_span=c(5, 9)
## freq_span=c(0, 25)
## amp_span=c(1e-06, 3e-05)
## cluster_span=c(4, 100)
## hhspec_image(hspec, time_span, freq_span, amp_span, cluster_span)


###################################################
### code chunk number 13: eemd_examples.Rnw:424-436 (eval = FALSE)
###################################################
## data(port_foster_event)
## ft=list()
## ft$xt=sig
## ft$dt=dt
## ft$nfft=4096
## ft$ns=30
## ft$nov=29
## 
## time_span=c(5, 10)
## freq_span=c(0, 25)
## amp_span=c(0.00001, 0.001)
## ftspec_image(ft, time_span, freq_span, amp_span)


