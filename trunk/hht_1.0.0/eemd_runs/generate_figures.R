#EEMD

w=1000 #PNG Pixel width
h=1000 #PNG Pixel height

#fig:nonlinstaemd
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

png(file="../hht/vignettes/nonlinstaemd.png", width=w, height=h)
plot_imfs(emd_result, time_span, imf_list, os, res)
dev.off()
print("Generated figure nonlinstaemd")

#fig:nonlinstarand
set.seed(628)
sig2=sin(2*pi*tt+0.5*sin(pi*tt))*exp(-(tt-50)^2/100)+sin(pi*tt/10)+rnorm(length(sig), 0, 0.01)
emd_result=sig2imf(sig2, dt, emd_config)
time_span=c(0, -1)
imf_list=c(1:emd_result$nimf)
os=TRUE
res=FALSE
png(file="../hht/vignettes/nonlinstarand.png", width=w, height=h)
plot_imfs(emd_result, time_span, imf_list, os, res)
dev.off()
print("Generated figure nonlinstarand")

#fig:modemix
library(hht)
set.seed(628)
dt=0.01
tt=seq_len(2000)*dt
sig=sin(2*pi*tt) + 2*sin(12*pi*tt)+rnorm(length(tt), 0, 0.4)
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
time_span=c(0, 20)
imf_list=1:5
os=TRUE
res=FALSE
png(file="../hht/vignettes/modemix.png", width=w, height=h)
plot_imfs(emd_result, time_span, imf_list, os, res)
dev.off()
print("Generated figure modemix")

#fig:modemixht
emd_result=hhtransform(emd_result)
max_freq=10
freq_step=0.1
hspec=hh_render(emd_result, max_freq, freq_step)

time_span=c(0, 20)
freq_span=c(0, 10)
amp_span=c(0.75, 3)
png(file="../hht/vignettes/modemixht.png", width=w, height=h)
hhspec_image(hspec, time_span, freq_span, amp_span)
dev.off()
print("Generated figure modemixht")

#fig:nonlinstasift
trials=100
nimf=11
noise_amp=0.1
trials_dir="example_1"
eemd_result=eemd_compile(trials_dir, trials, nimf)
emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=5
emd_config$stop_rule="type5"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
resift_rule="max_var"
resift_result=eemd_resift(eemd_result, emd_config, resift_rule)
time_span=c(0, 100)
imf_list=c(4, 5, 6, 8)
os=TRUE
res=TRUE
fit_line=TRUE
png(file="../hht/vignettes/nonlinstasift.png", width=w, height=h)
plot_imfs(resift_result, time_span, imf_list, os, res)
dev.off()
print("Generated figure nonlinstasift")

################################
##fig:stackedeemdimfs
#trials=100
#nimf=9
#trials_dir="example_2"
#eemd_result=eemd_compile(trials_dir, trials, nimf)
#
#emd_config=list()
#emd_config$max_sift=200
#emd_config$max_imf=100
#emd_config$tol=5
#emd_config$stop_rule="type5"
#emd_config$boundary="wave"
#emd_config$sm="none"
#emd_config$spar=NA
#emd_config$weight=20
#
#resift_rule="max_var"
#resift_result=eemd_resift(eemd_result, emd_config, resift_rule)
#time_span=c(0, 20)
#imf_list=1:5
#os=TRUE
#res=FALSE
#png(file="../hht/vignettes/stackedeemdimfs.png", width=w, height=h)
#plot_imfs(resift_result, time_span, imf_list, os, res)
#dev.off()
#print("generated fig:stackedeemdimfs")

##################################
##fig:stackedeemdht
#trials_dir="example_2"
#trials=100
#nimf=9
#eemd_result=eemd_compile(trials_dir, trials, nimf)
#
#max_freq=10
#freq_step=0.1
#hspec=hh_render(eemd_result, max_freq, freq_step)
#
#time_span=c(0, 20)
#freq_span=c(0, 10)
#amp_span=c(0.75, 3)
#cluster_span=c(8, 100)
#png(file="../hht/vignettes/stackedeemdht.png", width=w, height=h)
#hhspec_image(hspec, time_span, freq_span, amp_span, cluster_span)
#dev.off()
#print("Generated stackedeemdht")
##########################################
#fig:pferesift
#
#trials=100
#nimf=8
#noise_amp=6.4e-07
#trials_dir="example_3"
#eemd_result=eemd_compile(trials_dir, trials, nimf)
#list()
#emd_config$max_sift=200
#emd_config$max_imf=100
#emd_config$tol=5
#emd_config$stop_rule="type5"
#emd_config$boundary="wave"
#emd_config$sm="none"
#emd_config$spar=NA
#emd_config$weight=20
#resift_rule="max_var"
#resift_result=eemd_resift(eemd_result, emd_config, resift_rule)
#time_span=c(5, 10)
#imf_list=1:3
#os=TRUE
#res=FALSE
#fit_line=TRUE
#png(file="../hht/vignettes/pferesift.png", width=w, height=h)
#plot_imfs(resift_result, time_span, imf_list, os, res, fit_line)
#dev.off()
##print("Generated fig:pferesift")
#
#################################
##fig:pfeht
#trials_dir="example_3"
#trials=100
#nimf=11
#
#eemd_result=eemd_compile(trials_dir, trials, nimf)
#
#max_freq=30
#freq_step=0.05
#hspec=hh_render(eemd_result, max_freq, freq_step)
#
#time_span=c(5, 9)
#freq_span=c(0, 25)
#amp_span=c(1e-06, 3e-05)
#cluster_span=c(4, 100)
#png(file="../hht/vignettes/pfeht.png", width=w, height=h)
#hhspec_image(hspec, time_span, freq_span, amp_span, cluster_span)
#dev.off()
#print("Generated fig:pfeht")
#
##################################
##fig:pfeft
#data(port_foster_event)
#ft=list()
#ft$xt=sig
#ft$dt=dt
#ft$nfft=4096
#ft$ns=30
#ft$nov=29
#
#time_span=c(5, 10)
#freq_span=c(0, 25)
#amp_span=c(0.00001, 0.001)
#png(file="../hht/vignettes/pfeft.png", width=w, height=h)
#ftspec_image(ft, time_span, freq_span, amp_span)
#dev.off()
#print("Generated fig:pfeft")

################ INTERESTING SIGNALS

#fig:stacked1imfs

library(hht)
dt=0.01
tt=seq_len(10000)*dt
sig=0.5*sin(pi*tt/10)+sin(2*pi*tt)+0.5*sin(10*pi*tt)+0.5*sin(40*pi*tt)
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
emd_result=hhtransform(emd_result)

time_span=c(0,40)
imf_list=1:4
os=TRUE
res=FALSE
png(file="../hht/vignettes/stacked1imfs.png", width=w, height=h)
plot_imfs(emd_result, time_span, imf_list, os, res)
dev.off()
print("Generated fig:stacked1imfs")

#fig:stacked1farhspec

max_freq=25
freq_step=0.1
hspec=hh_render(emd_result, max_freq, freq_step)
time_span=c(0, -1)
freq_span=c(0, 25)
amp_span=c(0, 1)
png(file="../hht/vignettes/stacked1farhspec.png", width=w, height=h)
hhspec_image(hspec, time_span, freq_span, amp_span, cex=1)
dev.off()
print("Generated stacked1farhspec")

#fig:stacked1ft

ft=list()
ft$xt=sig
ft$dt=dt
ft$nfft=2048
ft$ns=1000
ft$nov=998
amp_span=c(100,200)
png(file="../hht/vignettes/stacked1ft.png", width=w, height=h)
ftspec_image(ft, time_span, freq_span, amp_span, cex=1)
dev.off()
print("generated stacked1ft")

#fig:stacked2imf

library(hht)
dt=0.05
tt=seq_len(10000)*dt
sig=sin(2*pi*tt)+0.25*sin(4*pi*tt)
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
time_span=c(20, 30)
imf_list=c(1)
os=TRUE
res=FALSE
fit_line=TRUE
png(file="../hht/vignettes/stacked2imf.png", width=w, height=h)
plot_imfs(emd_result, time_span, imf_list, os, res, fit_line)
dev.off()
print("Generated stacked2imf")

##fig:stacked2ht

emd_result=hhtransform(emd_result)
max_freq=4
freq_step=0.1
imf_list=c(1)
hspec=hh_render(emd_result, max_freq, freq_step, imf_list)
png(file="../hht/vignettes/stacked2ht.png", width=w, height=h)
time_span=c(20, 30)
freq_span=c(0, 3)
amp_span=c(0, -1)
hhspec_image(hspec, time_span, freq_span, amp_span, cex=1)
dev.off()

print("Generated stacked2ht")

#fig:stacked2ft

ft=list()
ft$xt=sig
ft$dt=dt
ft$nfft=4096
ft$ns=2000
ft$nov=1998
time_span=c(0, -1)
freq_span=c(0, 3)
amp_span=c(90, 200)
png(file="../hht/vignettes/stacked2ft.png", width=w, height=h)
ftspec_image(ft, time_span, freq_span, amp_span, cex=1)
dev.off()
print("Generated stacked2ft")

#fig:
