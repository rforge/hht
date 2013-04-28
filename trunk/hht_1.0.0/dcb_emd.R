#Most of this code is copied verbatim from the EMD package written by Donghoh Kim and Hee-Seok Oh
#I modified the code by adding the type 3 stop rule (S Stoppage Criterion)
#http://dasan.sejong.ac.kr/~dhkim/software_dcb_emd.html 

dcb_extractimf <- function(residue, tt=NULL, tol=sd(residue)*0.1^2, max.sift=20, 

                        stoprule="type1", boundary="periodic", sm="none", spar=NA, S=5, weight=20, check=FALSE) {     

    

    if (boundary == "none")

        minextrema <- 4

    else

        minextrema <- 2

            

    if (sm == "spline" & is.null(spar)) stop("Provide the smoothing parameter of spline smoothing.\n")

    ndata <- length(residue); ndatam1 <- ndata - 1 

    if(is.ts(residue)) 

        residue <- as.numeric(residue) 

    

    if(is.null(tt)) tt <- 1:length(residue)

    

    emin <- emax <- em <- h <- imf <- NULL

    n2data <- 2*ndata; tt2 <- 1:n2data

    n3data <- 3*ndata; n3datam1 <- n3data-1; tt3 <- 1:n3data



    input <- residue; rangext <- range(residue)

    

    j <- 1

    s_count=0 #Number of times zero crossings and extrema have stayed the same

    prev_excross=0 #How many times the previous sift crossed zero or hit an extrema

    repeat {

        tmp <- extrema(input, ndata, ndatam1)

        

        if (tmp$nextreme <= minextrema) break

        

        if(j == 1 || boundary == "wave") {

            minindex <- unique(c(t(tmp$minindex))); minn <- length(minindex)

            maxindex <- unique(c(t(tmp$maxindex))); maxn <- length(maxindex)          

  

            extminindex <- minindex; extmaxindex <- maxindex

            tmpwavefreq1 <- diff(sort(c(tt[1], tt[minindex[1]], tt[maxindex[1]])))

            tmpwavefreq2 <- diff(sort(c(tt[minindex[minn]], tt[maxindex[maxn]], tt[ndata])))

                                         

            if(input[1] <= input[minindex[1]] && input[1] <= input[maxindex[1]]) {

                extminindex <- c(1, extminindex); wavefreq1 <- 2 * tmpwavefreq1[1] 

            } else if(input[1] >= input[minindex[1]] && input[1] >= input[maxindex[1]]) {

                extmaxindex <- c(1, extmaxindex); wavefreq1 <- 2 * tmpwavefreq1[1] 

            } else if(input[1] >= (input[minindex[1]] + input[maxindex[1]])/2) {

                wavefreq1 <- tmpwavefreq1[2]  + max(tmpwavefreq1[2], 2*tmpwavefreq1[1])        

            } else {

                wavefreq1 <- tmpwavefreq1[2]  + max(tmpwavefreq1[2], round(1.5*tmpwavefreq1[1]))

            }                 

                  

            if(input[ndata] <= input[minindex[minn]] && input[ndata] <= input[maxindex[maxn]]) {

                extminindex <- c(extminindex, ndata); wavefreq2 <- 2 * tmpwavefreq2[2] 

            } else if(input[ndata] >= input[minindex[minn]] && input[ndata] >= input[maxindex[maxn]]) {

                extmaxindex <- c(extmaxindex, ndata); wavefreq2 <- 2 * tmpwavefreq2[2] 

            } else if(input[ndata] >= (input[minindex[minn]] + input[maxindex[maxn]])/2) {

                wavefreq2 <- tmpwavefreq2[1]  + max(tmpwavefreq2[1], 2*tmpwavefreq2[2])        

            } else {

                wavefreq2 <- tmpwavefreq2[1]  + max(tmpwavefreq2[1], round(1.5*tmpwavefreq2[2]))

            } 



            extminn <- length(extminindex); extmaxn <- length(extmaxindex)

            extttminindex <- c(tt[extminindex[1]] - 4:1 * wavefreq1, tt[extminindex], tt[extminindex[extminn]] + 1:4 * wavefreq2)

            extttmaxindex <- c(tt[extmaxindex[1]] - 4:1 * wavefreq1, tt[extmaxindex], tt[extmaxindex[extmaxn]] + 1:4 * wavefreq2)          

            

            if(sm == "none") {        

                f <- splinefun(extttminindex, c(rep(input[extminindex[1]], 4), input[extminindex], rep(input[extminindex[extminn]], 4)))

                emin <- cbind(emin, f(tt))

                f <- splinefun(extttmaxindex, c(rep(input[extmaxindex[1]], 4), input[extmaxindex], rep(input[extmaxindex[extmaxn]], 4)))

                emax <- cbind(emax, f(tt))

            } else if (sm == "spline") {                

                f <- sreg(extttminindex, c(rep(input[extminindex[1]], 4), input[extminindex], rep(input[extminindex[extminn]], 4)), lambda = spar)

                llambda <- f$lambda * weight 

                f <- sreg(extttminindex, c(rep(input[extminindex[1]], 4), input[extminindex], rep(input[extminindex[extminn]], 4)), lambda = llambda)

                llambda <- f$lambda                

                emin <- cbind(emin, predict(f, tt))

                f <- sreg(extttmaxindex, c(rep(input[extmaxindex[1]], 4), input[extmaxindex], rep(input[extmaxindex[extmaxn]], 4)), lambda = spar)

                ulambda <- f$lambda * weight

                f <- sreg(extttmaxindex, c(rep(input[extmaxindex[1]], 4), input[extmaxindex], rep(input[extmaxindex[extmaxn]], 4)), lambda = ulambda)

                ulambda <- f$lambda                

                emax <- cbind(emax, predict(f, tt)) 

            }

            

            em <- cbind(em, (emin[,j] + emax[,j]) / 2)

            

            if(check){

                plot(tt, input, type="l", col=3, xlab="", ylab="", main=paste("Boundary = ", boundary, sep="")) 

                points(tt[unique(c(t(tmp$minindex)))], input[unique(c(t(tmp$minindex)))], col=4)

                points(tt[unique(c(t(tmp$maxindex)))], input[unique(c(t(tmp$maxindex)))], col=2)

                lines(tt, emin[,j], col=4)

                lines(tt, emax[,j], col=2)

                lines(tt, em[,j]); locator(1)

            }

        } else {

        if(boundary == "none") {

            minindex <- c(1, unique(c(t(tmp$minindex))), ndata)

            maxindex <- c(1, unique(c(t(tmp$maxindex))), ndata)  

        

            fmin <- splinefun(tt[minindex], input[minindex])

            emin <- cbind(emin, fmin(tt))

        

            fmax <- splinefun(tt[maxindex], input[maxindex])

            emax <- cbind(emax, fmax(tt))

            

            em <- cbind(em, (emin[,j] + emax[,j]) / 2)    

                   

            if(check){

                plot(tt, input, type="l", col=3, xlab="", ylab="", main=paste("Boundary = ", boundary, sep=""))#,  ylim=rangext) 

                points(tt[unique(c(t(tmp$minindex)))], input[unique(c(t(tmp$minindex)))], col=4)

                points(tt[unique(c(t(tmp$maxindex)))], input[unique(c(t(tmp$maxindex)))], col=2)

                lines(tt, emin[,j], col=4)

                lines(tt, emax[,j], col=2)

                lines(tt, em[,j]); locator(1)

            }            

        }   

        if(boundary == "symmetric" || boundary == "periodic") {

            if(boundary == "symmetric") {

                inputext <- c(rev(input[-1]), input, rev(input)[-1])   

                ttext <- c(tt[1] - rev(cumsum(diff(tt))), tt, tt[ndata] + cumsum(rev(diff(tt))))        

                tmp <- extrema(inputext, n3data - 2, n3data - 3) 

            } else if (boundary == "periodic") {

                inputext <- c(input[-ndata], input, input[-1])   

                ttext <- c(rev(tt[1] - cumsum(rev(diff(tt)))), tt, tt[ndata] + cumsum(diff(tt)))        

                tmp <- extrema(inputext, n3data - 2, n3data - 3)            

            }           

            minindex <- unique(c(t(tmp$minindex)))

            #minindex <- minindex[minindex <= n2data + minindex[1]]

            maxindex <- unique(c(t(tmp$maxindex)))

            #maxindex <- maxindex[maxindex <= n2data + maxindex[1]]   

                        

            if(sm == "none") {     

                fmin <- splinefun(ttext[minindex], inputext[minindex])

                tmpmin <- fmin(ttext[ndata:(n2data-1)])

                

                fmax <- splinefun(ttext[maxindex], inputext[maxindex])

                tmpmax <- fmax(ttext[ndata:(n2data-1)])

            } else if (sm == "spline") { 

                fmin <- sreg(ttext[minindex], inputext[minindex], lambda = spar)

                llambda <- fmin$lambda * weight 

                fmin <- sreg(ttext[minindex], inputext[minindex], lambda = llambda)

                tmpmin <- predict(fmin, ttext[ndata:(n2data-1)])          

                

                fmax <- sreg(ttext[maxindex], inputext[maxindex], lambda = spar)

                ulambda <- fmax$lambda * weight                 

                fmax <- sreg(ttext[maxindex], inputext[maxindex], lambda = ulambda)

                tmpmax <- predict(fmax, ttext[ndata:(n2data-1)])                        

            } 

         

            emin <- cbind(emin, tmpmin[1:ndata])

            emax <- cbind(emax, tmpmax[1:ndata])

            

            tmpm <- (tmpmin+tmpmax)/2; em <- cbind(em, tmpm[1:ndata])



            if(check){

                plot(ttext[ndata:(n2data-1)], inputext[ndata:(n2data-1)], type="l", col=3, xlab="", ylab="", main=paste("Boundary = ", boundary, sep=""))

                points(ttext[minindex[minindex > ndata & minindex < n2data]], inputext[minindex[minindex > ndata & minindex < n2data]], col=4)

                points(ttext[maxindex[maxindex > ndata & maxindex < n2data]], inputext[maxindex[maxindex > ndata & maxindex < n2data]], col=2)

                lines(ttext[ndata:(n2data-1)], tmpmin, col=4)

                lines(ttext[ndata:(n2data-1)], tmpmax, col=2)

                lines(ttext[ndata:(n2data-1)], tmpm); locator(1)

            } 

        } else if(boundary == "evenodd") {

            inputeven <- c(input, rev(input), input) 

            ttext <- c(tt, tt[ndata]+tt[ndata]-tt[ndatam1], tt[ndata]+tt[ndata]-tt[ndatam1] + cumsum(rev(diff(tt))))

          

            tmp <- extrema(inputeven, n3data, n3datam1) 



            minindex <- unique(c(t(tmp$minindex)))

            minindex <- minindex[minindex <= n2data + minindex[1]]

            maxindex <- unique(c(t(tmp$maxindex)))

            maxindex <- maxindex[maxindex <= n2data + maxindex[1]]

            

            if(sm == "none") { 

                f <- splinefun(ttext[minindex], inputeven[minindex])

                emineven <- f(ttext[1:ndata]) 

                #lines(tt2, emineven, col=4)



                f <- splinefun(ttext[maxindex], inputeven[maxindex])

                emaxeven <- f(ttext[1:ndata])

                #lines(tt2, emaxeven, col=2)

            } else if (sm == "spline") { 

                fmin <- sreg(ttext[minindex], inputeven[minindex], lambda = spar)

                llambda <- fmin$lambda * weight 

                fmin <- sreg(ttext[minindex], inputeven[minindex], lambda= llambda)

                emineven <- predict(fmin, ttext[1:ndata])          

                

                fmax <- sreg(ttext[maxindex], inputext[maxindex], lambda = spar)

                ulambda <- fmax$lambda * weight                 

                fmax <- sreg(ttext[maxindex], inputext[maxindex], lambda = ulambda)

                emaxeven <- predict(fmax, ttext[1:ndata])                        

            }             

            

            inputodd <- c(input, -rev(input), input)

            tmp <- extrema(inputodd, n3data, n3datam1) 



            minindex <- unique(c(t(tmp$minindex)))

            minindex <- minindex[minindex <= n2data + minindex[1]]

            maxindex <- unique(c(t(tmp$maxindex)))

            maxindex <- maxindex[maxindex <= n2data + maxindex[1]]            



            if(sm == "none") { 

                f <- splinefun(ttext[minindex], inputodd[minindex])

                eminodd <- f(ttext[1:ndata]) 

                #lines(tt2, eminodd, col=4)



                f <- splinefun(ttext[maxindex], inputodd[maxindex])

                emaxodd <- f(ttext[1:ndata])

                #lines(tt2, emaxodd, col=2)

            } else if (sm == "spline") { 

                fmin <- sreg(ttext[minindex], inputodd[minindex], lambda = spar)

                llambda <- fmin$lambda * weight 

                fmin <- sreg(ttext[minindex], inputodd[minindex], lambda = llambda)

                eminodd <- predict(fmin, ttext[1:ndata])          

                

                fmax <- sreg(ttext[maxindex], inputext[maxindex], lambda = spar)

                ulambda <- fmax$lambda * weight                 

                fmax <- sreg(ttext[maxindex], inputext[maxindex], lambda = ulambda)

                emaxodd <- predict(fmax, ttext[1:ndata])                        

            }  

  

            emin <- cbind(emin, (emineven+eminodd)/2)

            emax <- cbind(emax, (emaxeven+emaxodd)/2)

            em <- cbind(em, (emin[,j]+emax[,j])/2)

            

            if(check) {

                plot(tt, (inputeven+inputodd)/2[1:ndata], type="l", col=3, ylab="")

                lines(tt, emin, col=4)

                lines(tt, emax, col=2)

                lines(tt, (emin+emax)/2); locator(1)

            }  

        }

        }

        

        h <- cbind(h, input - em[,j])    



        if (stoprule == "type1" && (all(abs(em[,j]) < tol))) { 

            imf <- h[,j]

            residue <- residue - imf

            break

        }

                        

        if (stoprule == "type2" && j >= 2) {

            if (sum((h[2:ndatam1, j-1]-h[2:ndatam1, j])^2/h[2:ndatam1, j-1]^2) < tol) { 

                imf <- h[,j]

                residue <- residue - imf

                break

                }

        }

	if(stoprule=="type3" && j >= 2) #If the number of extrema and zero crossings do not change with successive siftings

        {

		if(tmp$nextreme+tmp$ncross==prev_excross) #If the extrema/zero crossings have not changed then

		{

			s_count=s_count+1

			if (s_count>=S) #If the IMF has not changed in S iterations, break and return IMF

			{

				imf=h[,j]

				residue <- residue - imf

				break

			}	

		}

		prev_excross=tmp$nextreme+tmp$ncross

	}

        

	if (j+1>max.sift)

	{

		imf=h[,j]

		residue=residue - imf

		print("Exceeded maximum sifts!")

		break

	}

        input <- h[,j]

       

        j <- j+1

    }

    

    if(check) list(emin=emin, emax=emax, em=em, h=h, imf=imf, residue=residue, niter=j) else   

    list(imf=imf, residue=residue, niter=j)

}





dcb_emd <- function(xt, tt=NULL, tol=sd(xt)*0.1^2, max.sift=20, stoprule="type1", boundary="periodic", 

                smlevels=c(1), sm="none", spar=NA, weight=20, 

                check=FALSE, max.imf=10, plot.imf=TRUE, interm=NULL,S=5) {



    sifts=NULL #Count how many sifts each IMF has


    if(is.ts(xt))

        xt <- as.numeric(xt) 



    if(is.null(tt)) tt <- 1:length(xt)

        

    if(is.null(interm) || all(interm <= 0)) intermtest <- FALSE else intermtest <- TRUE

    ndata <- length(xt); ndatam1 <- ndata - 1

    residue <- xt; rangext <- range(residue)

    imf <- NULL

    

    j <- 1

    

    #firstimf <- dcb_extractimf(residue, tt, tol, max.sift, 

    #                        stoprule=stoprule, boundary=boundary, check=check)$imf

    #if(!is.null(firstimf)) rangeimf <- range(firstimf)

    

    repeat {

        if (j > max.imf) break

        if ((any(j == smlevels) || smlevels == "all") & sm == "spline")

        {

            tmp <- dcb_extractimf(residue, tt, tol, max.sift, 

                            stoprule=stoprule, boundary=boundary, sm=sm, spar=spar, check=check) 

	    sifts=append(sifts,tmp$niter)

	}

	else

	{

            tmp <- dcb_extractimf(residue, tt, tol, max.sift, 

                            stoprule=stoprule, boundary=boundary, check=check)

            sifts=append(sifts,tmp$niter)

	}

#        if(j == 1 && !is.null(tmp$imf)) rangeimf <- range(tmp$imf)         

#        if (tmpstop$nextreme == 0 || tmpstop$ncross == 0 || tmpstop$nextreme != tmpstop$ncross ||

#            tmpstop$nextreme != (tmpstop$ncross+1) || j >= max.imf) {

#            break

#        }

        if (is.null(tmp$imf)) {

            break

        }

        

        if(plot.imf) {

            plot(tt, residue, type="l", xlab="", ylab="", #ylim=rangext,

                main=paste(j-1, "-th residue=", j, "-th imf+", j, "-th residue", sep="")); abline(h=0)

        }     



        if(intermtest && length(interm) >= j && interm[j] > 0) {

            

            tmpimf <- tmp$imf

            tmpresidue <- tmp$residue



            tmpinterm <- extrema(tmpimf, ndata, ndatam1)

            tmpncross <- tmpinterm$ncross

            zerocross <- as.numeric(round(apply(tmpinterm$cross, 1, mean)))



            if(abs(tt[zerocross[3]] - tt[zerocross[1]]) > interm[j]) {

                tmpresidue[1:zerocross[3]] <- tmpresidue[1:zerocross[3]] + tmpimf[1:zerocross[3]]                     

                tmpimf[1:zerocross[3]] <- 0            

            }

            

            for (k in seq(3, tmpncross-3, by=2))         

                if(abs(tt[zerocross[k+2]] - tt[zerocross[k]]) > interm[j]) {

                    tmpresidue[zerocross[k]:zerocross[k+2]] <- tmpresidue[zerocross[k]:zerocross[k+2]] + 

                                tmpimf[zerocross[k]:zerocross[k+2]]            

                    tmpimf[zerocross[k]:zerocross[k+2]] <- 0

                }

            

            if(!(tmpncross %% 2)) {               

                if(abs(tt[zerocross[tmpncross]] -tt[zerocross[tmpncross-1]]) > interm[j]/2) {       

                    tmpresidue[zerocross[tmpncross-1]:ndata] <- 

                            tmpresidue[zerocross[tmpncross-1]:ndata] + tmpimf[zerocross[tmpncross-1]:ndata]                 

                    tmpimf[zerocross[tmpncross-1]:ndata] <- 0            

                } 

            } else {

                if(abs(tt[zerocross[tmpncross]] - tt[zerocross[tmpncross-2]]) > interm[j])       

                    tmpresidue[zerocross[tmpncross-2]:ndata] <- 

                            tmpresidue[zerocross[tmpncross-2]:ndata] + tmpimf[zerocross[tmpncross-2]:ndata]                 

                    tmpimf[zerocross[tmpncross-2]:ndata] <- 0              

            }

            

            tmp$imf <- tmpimf

            tmp$residue <- tmpresidue            

        }

        

        imf <- cbind(imf, tmp$imf)     

        residue <- tmp$residue     

             

        if(plot.imf) {

            plot(tt, imf[,j], type="l", xlab="", ylab="", #ylim=rangeimf, 

                main=paste(j, "-th imf", sep="")); abline(h=0)

            plot(tt, residue, type="l", xlab="", ylab="", #ylim=rangext, 

                main=paste(j, "-th residue", sep="")); abline(h=0); locator(1)

        }

        print(paste("Extracted IMF",as.character(j)))

        j <- j+1

    }

    list(imf=imf, residue=residue, nimf=j-1,sifts=sifts)

}  


