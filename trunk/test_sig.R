source("hht/R/hht.R")
library(fields)

trials_dir = "/home/bowman/research/hht/figures/srl_computational/type2_eemd"
trials = 100
nimf = 10
eemd_result = eemd_compile(trials_dir, trials, nimf)

dt = 0.01
dfreq = 0.01
freq_span = NULL
scaling = "log"
verbose = TRUE

hspec = hh_render(eemd_result, dt, dfreq, freq_span = freq_span, scaling = scaling, verbose = verbose)

hhspec_image(hspec, time_span = c(4, 10), freq_span = c(0, 10), clusterspec = FALSE, img_y_lab = "Log Frequency", pretty = TRUE, scaling = "log", amp_span = c(1e-6, 1e-4), cluster_span =c(5, 100))
