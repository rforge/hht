
library(EMD)
############  EEMD function

eemd <- function(
xt, tt=NULL, sdn=0.2, n=100, 
max.sift=20, stoprule="type1", boundary="periodic", max.imf=10)
{

    ndata <- length(xt); noiselevel <- sd(xt)*sdn

    imf <- matrix(0, ndata, max.imf)
    residue <- rep(0, ndata)
    
    for (i in 1:n) {
        yt <- xt + rnorm(ndata, 0, noiselevel)
        tmp <- emd(xt=yt, tt=tt, tol=sd(yt)*0.1^2, max.sift=max.sift, stoprule=stoprule,
                   boundary=boundary, sm="none", smlevels=c(1), spar=NULL, alpha=NULL, 
                   check=FALSE, max.imf=max.imf, plot.imf=FALSE, interm=NULL, weight=NULL)
        for (j in 1:tmp$nimf)
            imf[, j] <- imf[, j] + tmp$imf[, j]
           
        residue <- residue + tmp$residue
    }
    
   list(imf=imf/n, residue=residue/n)
}



### EMD and EEMD for noisy signal
ndata <- 2048                                                                      
                                     
tt3 <- seq(0, 9, length=ndata)                                                     
xt3 <- sin(pi * tt3) + sin(2* pi * tt3) + sin(6 * pi * tt3)  + 0.5 * tt3    
set.seed(1)
xt4 <- xt3 + rnorm(ndata, 0, sd(xt3)/5)

emdresult <- emd(xt4, tt3, boundary = "wave", max.imf = 4)                                     

set.seed(77)
eemdresult <- eemd(xt4, tt3, sdn=0.3, n=100,  boundary="wave", max.imf=7)

par(mfcol=c(9,1), mar=c(2,2,2,1))                                                 
rangext <- range(xt3); rangeimf <- rangext - mean(rangext)

plot(tt3, xt4, xlab="", ylab="", main="Ensemble EMD", ylim=rangext, type="l")                      
plot(tt3, eemdresult$imf[,1], xlab="", ylab="", main="imf 1", ylim=rangeimf,  type="l"); abline(h=0, lty=2) 
plot(tt3, eemdresult$imf[,2], xlab="", ylab="", main="imf 2", ylim=rangeimf,  type="l"); abline(h=0, lty=2) 
plot(tt3, eemdresult$imf[,3], xlab="", ylab="", main="imf 3", ylim=rangeimf,  type="l"); abline(h=0, lty=2)
plot(tt3, eemdresult$imf[,4], xlab="", ylab="", main="imf 4", ylim=rangeimf,  type="l"); abline(h=0, lty=2)
plot(tt3, eemdresult$imf[,5], xlab="", ylab="", main="imf 5", ylim=rangeimf,  type="l"); abline(h=0, lty=2)
plot(tt3, eemdresult$imf[,6], xlab="", ylab="", main="imf 6", ylim=rangeimf,  type="l"); abline(h=0, lty=2)
plot(tt3, eemdresult$imf[,7], xlab="", ylab="", main="imf 6", ylim=rangeimf,  type="l"); abline(h=0, lty=2)
plot(tt3, eemdresult$residue, xlab="", ylab="", main="residue", ylim=rangext, type="l") 
