source("hht/R/hht.R")
library(fields)

trials.dir = "/home/bowman/research/hht/figures/srl.computational/type2.EEMD"
trials = 100
nimf = 10
EEMD.result = EEMDCompile(trials.dir, trials, nimf)

dt = 0.01
dfreq = 0.01
freq.span = NULL
scaling = "none"
verbose = TRUE

hgram = HHRender(EEMD.result, dt, dfreq, freq.span = freq.span, scaling = scaling, verbose = verbose)

HHGramImage(hgram)
