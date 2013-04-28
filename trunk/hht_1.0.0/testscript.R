#Test EMD
source("hht.R")
source("dcb_emd.R")
require(EMD)
require(RSEIS)

dt=0.01
tt=seq_len(10000)*dt
sig=sin(2*pi*tt)+0.5*sin((pi/10)*tt)+0.25*sin(10*pi*tt)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type3"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_config$S=5

emd_result=sig2imf(sig, dt, emd_config)
emd_result=hhtransform(emd_result)
max_freq=6
freq_step=0.01
hspec=hh_render(emd_result, max_freq, freq_step)

time_span=c(0, -1)
freq_span=c(0, 6)
amp_span=c(0, -1)

#hhspec_image(hspec, time_span, freq_span, amp_span)

#Test EEMD
dt=0.01
tt=seq_len(10000)*dt
sig=sin(2*pi*tt)+0.5*sin((pi/10)*tt)+0.25*sin(10*pi*tt)

emd_config=list()
emd_config$max_sift=200
emd_config$max_imf=100
emd_config$tol=0.2
emd_config$stop_rule="type3"
emd_config$boundary="wave"
emd_config$sm="none"
emd_config$spar=NA
emd_config$weight=20
emd_config$S=5

trials=10
nimf=10
noise_amp=0.1
trials_dir="test"

eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)

eemd_result=eemd_compile(trials_dir, trials, nimf)

time_span=c(0, -1)
imf_list=1:9
os=TRUE
res=TRUE
plot_imfs(eemd_result, time_span, imf_list, os, res)

max_freq=6
freq_step=0.01
hspec=hh_render(eemd_result, max_freq, freq_step)

save(hspec, file="hspec.RDATA")
time_span=c(0, -1)
freq_span=c(0, 6)
amp_span=c(0, 1)

hhspec_image(hspec, time_span, freq_span, amp_span)

#Test resifting

resift_rule="max_var"

resift_result=eemd_resift(eemd_result, emd_config, resift_rule)
resift_result=hhtransform(resift_result)
hspec=hh_render(resift_result, max_freq, freq_step)
plot_imfs(resift_result, time_span, imf_list, os, res)
hhspec_image(hspec, time_span, freq_span, amp_spanm
