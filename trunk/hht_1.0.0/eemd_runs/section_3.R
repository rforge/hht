library(hht)
data(port_foster_event)
set.seed(628)
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
noise_amp=6.4e-07
trials_dir="example_3"

eemd(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir)
eemd_result=eemd_compile(trials_dir, trials, nimf)
