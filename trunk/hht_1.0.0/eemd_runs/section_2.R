library(hht)
set.seed(628)
dt=0.01
tt=seq_len(2000)*dt
sig=2*sin(2*pi*tt) + 0.5*sin(12*pi*tt)+rnorm(length(tt), 0, 0.4)
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
nimf=8
noise_amp=0.4
trials_dir="example_2"

eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)
eemd_result=eemd_compile(trials_dir, trials, nimf)
