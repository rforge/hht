library(hht)
source("hht/R/empirical_mode_decomposition.R")

set.seed(628)

noise.amp <- 0.1
tt <- seq(0, 10, by = 0.01)
sig <- 2 * sin(2 * pi * tt) + sin(10 * pi * tt) + noise.amp * rnorm(length(tt))

plot(tt, sig, type = "l")

emd.result <- Sig2IMF(sig, tt)

PlotIMFs(emd.result)

trials <- 100

verbose = TRUE
spectral.method = "arctan"
diff.lag = 1
 tol = 5
 max.sift = 200
stop.rule = "type5"
boundary = "wave"
sm = "none"
smlevels = c(1)
spar = NULL
max.imf = 100
interm = NULL
noise.type = "gaussian"
noise.array = NULL



ceemd.result <- CEEMD(sig, tt, noise.amp, 100)
dev.new()
PlotIMFs(ceemd.result)

#Try with Deception Island

data(PortFosterEvent)
noise.amp <- 6.4e-07
trials <- 100

di.emd <- Sig2IMF(sig, tt)
PlotIMFs(di.emd)

di.ceemd <- CEEMD(sig, tt, noise.amp, trials)
dev.new()
PlotIMFs(di.ceemd)
