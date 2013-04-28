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


trials=100
nimf=11
noise_amp=0.1
trials_dir="example_1"

eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)
eemd_result=eemd_compile(trials_dir, trials, nimf)
