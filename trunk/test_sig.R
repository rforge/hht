source("hht/R/hht.R")
library(fields)

trials_dir = "/home/bowman/research/hht/figures/srl_computational/type2_eemd"
trials = 100
nimf = 10
eemd_result = eemd_compile(trials_dir, trials, nimf)

dt = 0.01
dfreq = 0.01
freq_span = NULL
scaling = "none"
verbose = TRUE

hgram = hh_render(eemd_result, dt, dfreq, freq_span = freq_span, scaling = scaling, verbose = verbose)

hhgram_image(hgram)
