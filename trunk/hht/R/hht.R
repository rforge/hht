combine_trials <- function(in_dirs, out_dir, copy = TRUE)
{
   #Moves trial files from different directories, numbers them sequentially, and puts them in the specified directory.
   #This is important because the function eemd_compile expects them to be numbered consecutively from 1, and will crash if this is not the case.
   #INPUTS
   #    IN_DIRS is a vector of paths to directories containing trial files.  Trial files will be found recursively (so trial files in subdirectories will be discovered and moved).
   #    OUT_DIR is a directory to put the combined trial set into.  If OUT_DIR does not exist, it will be created
   #    COPY asks if you want to copy files into OUT_DIR (TRUE) or move them from IN_DIRS to OUT_DIR (false)

   trial_file_pattern = "TRIAL_\\d+\\.?(RData|RDATA)$" 
   e1=simpleError(paste("Directory", out_dir, "is not empty!"))
   if(!file.exists(out_dir))
   {
       dir.create(out_dir, recursive = TRUE)
       cat("Created directory:", out_dir, "\n")
   }

   if(length(list.files(out_dir))>0)
   {
       stop(e1)
   }

   c = 1
   for(d in in_dirs)
   {
       if(!file.exists(d)) 
       {
           warning(cat("Trials directory", d, "does not exist and will be skipped."))
       }
       else
       {
            trial_files = list.files(d, pattern = trial_file_pattern, recursive = TRUE, full.names = TRUE)
            if(length(trial_files) == 0)
            {
                warning(cat("No EMD trial files found in directory", d))
            }
            else
            {
               for (trial_file in trial_files)
               {
                   if(copy)
                   {
                       res = file.copy(trial_file, paste(out_dir, "/", "TRIAL_",sprintf("%05i",c),".RData", sep = ""))
                       if(!res)
                       {
                           warning(cat("Failed to copy", trial_file))
                       }
                   }
                   else
                   {
                       file.rename(trial_file, cat(out_dir, "/", "TRIAL_",sprintf("%05i",c),".RData", sep = ""))
                   }

                   c = c + 1
               }
           }
       }
   }
}

cosine_taper <- function(x, taper = 0.1)
{
    #Copied verbatim from the RSEIS package ver. 3.0-6 by Jonathan M. Lees
    #INPUTS
    #    X is the signal to taper
    #    TAPER is the fraction of the signal to apply tapering to (greater than 0, less than 0.5) 
    #OUTPUTS
    #    X is the tapered signal
    if (any(taper < 0) || any(taper > 0.5)) 
        stop("'taper' must be between 0 and 0.5")
    a <- attributes(x)
    x <- as.matrix(x)
    nc <- ncol(x)
    if (length(taper) == 1) 
        taper <- rep(taper, nc)
    else if (length(taper) != nc) 
        stop("length of 'taper' must be 1 or equal the number of columns of 'x'")
    nr <- nrow(x)
    for (i in 1:nc) {
        m <- floor(nr * taper[i])
        if (m == 0) 
            next
        w <- 0.5 * (1 - cos(pi * seq.int(1, 2 * m - 1, by = 2)/(2 * 
            m)))
        x[, i] <- c(w, rep(1, nr - 2 * m), rev(w)) * x[, i]
    }
    attributes(x) <- a
    x
}

eemd <-function(sig, dt, trials, nimf, noise_amp, emd_config, trials_dir=NULL)
{
	#Performs the Ensemble Empirical Mode Decomposition as described in Huang and Wu (2008) A Review on the Hilbert Huang Transform Method and its Applications to Geophysical Studies
	#It runs EMD on a given signal for N=TRIALS
 	#Each IMF set is saved to disk in TRIALS_DIR, which is created if it does not exist already.
	#Finally the EEMD function averages IMFs from all the trials together to produce an ensemble average and saves it in TRIALS_DIR 
	#INPUTS
	#	SIG is the time series to be analyzed
	#	DT is the sample rate
	#	TRIALS is the number of EMD analyses to be run
	#	NIMF is the number of IMFs to save; IMFs past this number will not be averaged
        #       NOISE_AMP determines what amplitude to make the white noise
	#	EMD_CONFIG determines how the EMD algorithm operates.
	#		EMD_CONFIG$MAX_SIFT how many times the IMFs can be sifted
	#		EMD_CONFIG$MAX_IMF maximum number of IMFs that can be generated
	#    	   	EMD_CONFIG$TOL Sifting stop criterion.
        #		EMD_CONFIG$STOP_RULE Make sure to read section on stop rules and make sure you understand what they imply!
        #		EMD_CONFIG$BOUNDARY How the start and stop of the time series are handled duing the spline process.
        #		EMD_CONFIG$SM Spline smoothing
        #		EMD_CONFIG$SPAR Smoothing parameter (only needed if sm is not none)
        #		EMD_CONFIG$WEIGHT Weight if "sm" is "spline"
	#	TRIALS_DIR is the location to store files generated during EEMD trials, if NULL then this program creates
	#	a directory called "trials" in the current directory
	#OUTPUTS are saved to TRIALS_DIR in variable EMD_RESULT
	
	if(is.null(trials_dir))
	{
		trials_dir="trials"
	}

	if(!file.exists(trials_dir))
	{
		dir.create(trials_dir, recursive=TRUE)
		print(paste("Created trial directory:",trials_dir))
	}

	averaged_imfs=array(0,nimf*length(sig),dim=c(length(sig),nimf))
	averaged_noise=array(0,length(sig),dim=c(length(sig),1))
	averaged_residue=array(0,length(sig),dim=c(length(sig),1))
	for (j in 1:trials)
	{
	
		noise=runif(length(sig),min=noise_amp*-1, max=noise_amp)
		tmpsig=sig+noise
                tt=seq_len(length(sig))*dt
		emd_result=sig2imf(tmpsig,dt, emd_config)
		emd_result$noise=noise
		emd_result$original_signal=tmpsig-noise
		save(emd_result, file=paste(trials_dir, "/", "TRIAL_",sprintf("%05i",j),".RData",sep=""))
		if (emd_result$nimf<nimf)
		{
			trial_nimf=emd_result$nimf
		}
		else
		{
			trial_nimf=nimf
		}
		averaged_imfs[,1:trial_nimf]=averaged_imfs[,1:trial_nimf]+emd_result$imf[,1:trial_nimf]
		averaged_noise=averaged_noise+noise
		averaged_residue=averaged_residue+emd_result$residue
		print(paste("TRIAL",as.character(j),"OF",as.character(trials),"COMPLETE"))
	}
}

eemd_compile<-function(trials_dir, trials, nimf)
{
	#Averages trials together to produce a set of ensemble IMFs.
	#Produces the Hilbert spectrogram of each trial and puts it into a structure with all the other trials.
	#This is used later to generate an ensemble Hilbert spectrogram of the entire EEMD run.
	#INPUTS
	#	TRIALS_DIR is the location where the trial files produced by EEMD are stored
	#	TRIALS is the number of trials to average together
	#	NIMF is the number of IMFs to build
	#OUTPUTS
	#	EEMD_RESULT containes the ensemble IMFs

        trial_file_pattern = "TRIAL_\\d+\\.?(RData|RDATA)$"
	emd_result=NULL
	if(length(list.files(trials_dir, pattern = trial_file_pattern))==0)
	{
		stop(paste("No EMD trial files found in directory", trials_dir))
	}
			
	counter=1
	sind=1
        for(file_name in list.files(trials_dir, pattern=trial_file_pattern))
        {
		load(paste(trials_dir,"/",file_name,sep=""))
		
		
                if(counter==1)
                {
			siglen=length(emd_result$original_signal)
		        averaged_imfs=array(0,nimf*siglen,dim=c(siglen,nimf))
        		averaged_noise=array(0,siglen,dim=c(siglen,1))
        		averaged_residue=array(0,siglen,dim=c(siglen,1))
			hinstfreq=array(rep(0,length(emd_result$original_signal)*nimf*trials),dim=c(length(emd_result$original_signal),nimf,trials))
			hamp=array(rep(0,length(emd_result$original_signal)*nimf*trials),dim=c(length(emd_result$original_signal),nimf,trials))
                }
                if(emd_result$nimf>=nimf)
                {
                        imf_ind=nimf
                }
                else
                {
                        imf_ind=emd_result$nimf
                }
		averaged_imfs[,1:imf_ind]=averaged_imfs[,1:imf_ind]+emd_result$imf[,1:imf_ind]
		averaged_noise=averaged_noise+emd_result$noise
		averaged_residue=averaged_residue+emd_result$residue

                hinstfreq[,1:imf_ind,counter]=emd_result$hinstfreq[,1:imf_ind]
                hamp[,1:imf_ind,counter]=emd_result$hamp[,1:imf_ind]
                sind=sind+imf_ind

                counter=counter+1
                if(counter>trials)
                {
                        break
                }
        }
	counter=counter-1
        real_imfs = which(apply(array(as.logical(averaged_imfs), dim = dim(averaged_imfs)), 2, any) == TRUE) #Find out which IMFs have data in them
        if(length(real_imfs < nimf))
        {
            warning("The number of requested IMFs is greater than the maximum number of IMFs produced during individual EEMD trials.  Only IMFs with data will be recorded.")
        }

	averaged_imfs=averaged_imfs[, real_imfs] / counter
        hinstfreq = hinstfreq[, real_imfs, ]
        hamp = hamp[, real_imfs, ]
	averaged_noise=averaged_noise/counter
	averaged_residue=averaged_residue/counter
	
        if(counter<trials)
        {
                warning("Number of trials requested was greater than the number of trials found in the trials directory")
        }

        
	eemd_result=c()
        eemd_result$nimf=length(real_imfs)
        eemd_result$dt=emd_result$dt
        eemd_result$original_signal=emd_result$original_signal
	eemd_result$averaged_imfs=averaged_imfs
	eemd_result$averaged_noise=averaged_noise
	eemd_result$averaged_residue=averaged_residue
        eemd_result$hinstfreq=hinstfreq
        eemd_result$hamp=hamp
	eemd_result$trials=trials 
        invisible(eemd_result)
}

eemd_resift <- function(eemd_result, emd_config, resift_rule)
{
	#Resifts averaged IMFs generated by EEMD to generate valid IMFs for Hilbert Transform
	#INPUTS
	#	EEMD_RESULT contains the averaged IMF set generated from EEMD.
	#	EMD_PARAMS controls how the EMD of the averaged IMFs is handled.
        #         EMD_CONFIG$MAX_SIFT how many times the IMFs can be sifted
        #         EMD_CONFIG$MAX_IMF maximum number of IMFs that can be generated
        #         EMD_CONFIG$TOL Sifting stop criterion.
        #         EMD_CONFIG$STOP_RULE Make sure to read section on stop rules and make sure you understand what they imply!
        #         EMD_CONFIG$BOUNDARY How the start and stop of the time series are handled duing the spline process.
        #         EMD_CONFIG$SM Spline smoothing
        #         EMD_CONFIG$SPAR Smoothing parameter (only needed if sm is not none)
        #         EMD_CONFIG$WEIGHT Weight if "sm" is "spline"
	#	RESIFT_RULE determines how the resifting occurs
	#	If resift_rule is numeric, get the nth IMF (so if resift_rule is 2, get the 2nd resifted IMF)
	#	If resift_rule is "last", return the last IMF.
	#	If resift_rule is "max_var", return the IMF with the most variance
	#	If resift_rule is "all", get all the IMFs returned by rerunning EMD on the averaged IMFs made by EEMD.
	#	This will likely be quite large.
	#OUTPUTS
	#	EEMD_RESULT$IMF a set of IMFs that are generated from the EEMD imfs.

	resift_result=eemd_result
	resift_result$imf=c()
	resift_result$averaged_imfs=NULL

	if(!is.numeric(resift_rule) & !resift_rule %in% c("last", "max_var", "all"))
        {
               e=simpleError(paste("Did not recognize resift_rule:", resift_rule))
               stop(e)
        }		
	for(k in seq_len(dim(eemd_result$averaged_imfs)[2]))
	{
		if(sum(eemd_result$averaged_imfs[,k]==0)!=length(eemd_result$averaged_imfs[,k]))
		{
			emd_result=sig2imf(eemd_result$averaged_imfs[,k], eemd_result$dt, emd_config)
			if(is.numeric(resift_rule))
			{
				if(emd_result$nimf>=resift_rule)
				{
					resift_result$imf=cbind(resift_result$imf, emd_result$imf[,resift_rule])
				}
				else
				{
					resift_result$imf=cbind(resift_result$imf, NA)
				}
			}
			else
			{
				if(resift_rule=="last")
				{
					resift_result$imf=cbind(resift_result$imf, emd_result$imf[,emd_result$nimf])
				}
				
				if(resift_rule=="max_var")
				{
					var_list=c()
					for(j in seq_len(emd_result$nimf))
					{
						var_list=cbind(var_list, var(emd_result$imf[,j]))
					}
					
					resift_result$imf=cbind(resift_result$imf, emd_result$imf[,var_list==max(var_list)])
				}

				if(resift_rule=="all")
				{
					resift_result$imf=cbind(resift_result$imf, emd_result$imf)
				}
			}	
		}
	}
	
	resift_result$resift_emd_config=emd_config
	resift_result$resift_rule=resift_rule
	resift_result$nimf=dim(resift_result$imf)[2]
	resift_result=hhtransform(resift_result)
	invisible(resift_result)
}

evolutive_fft <- function(sig, dt, ft, freq_span, taper = 0.05)
{
    #Calculates the evolutive Fourier spectrogram for use in FTSPEC_IMAGE
    #This code is modified the evolfft function in the RSEIS package.
    #INPUTS
    #    SIG is the signal to analyze
    #    DT is the sample rate
    #    FT is the fourier transform input parameters
    #       FT$NFFT is the fft length
    #       FT$NS is the number of samples in a window
    #       FT$NOV is the number of samples to overlap
    #    FREQ_SPAN is the frequency range to return
    #    TAPER is the percent taper applied to each window
    #OUTPUTS
    #    RET is the spectrogram image

    NT = length(sig);
    Nfft=ft$nfft
    Ns = ft$ns
    Nov = ft$nov
    nyquistf = 1/(2*dt);
    if(Nov<1)
      {
        Nov = floor(Ns - 0.1*Ns);
      }

    Ns = floor(Ns)

    if(Ns>NT)
      {
        emsg = c("ERROR: illegal call to evolutive_fft.",
        "Number of samples in trace must be greater than the number of sample in the moving window")
        cat(emsg, sep="\n")
        return(NULL)
      }

    kcol =floor( (NT-floor(Nov) )/(Ns-floor(Nov)))
    min1 = Nfft%%2;
    if(min1 == 0)
      {
        ## /* even */
        krow = (Nfft/2);
      } else {
        ##  /*  odd */
        krow = (Nfft+1)/2;
      }

    if(krow - Ns < 0) 
    {
        print("Error in evolfft: The number of rows in the spectrogram matrix is less than the number of samples in the window.  Increase Nfft and try again.")
        return(NULL)
    }
    skiplen = Ns - Nov;

    df = 1.0/(Nfft*dt);
    numfreqs=krow;
    if(kcol<1)
      {
        cat(paste(sep=' ', "error in evolfft kcol=", kcol, "krow=", krow, "NT", NT, "Ns", Ns, "Nov", Nov), sep="\n")
        return(NULL)
      }

    DMAT = matrix(rep(0,krow*kcol), ncol=kcol, nrow=krow)

    m = 1:(kcol)
    ibeg=((m-1)*skiplen)+1;
     iend = floor(ibeg+Ns-1)

    for( i in m)
      {
        tem = sig[ibeg[i]:iend[i]]
        tem = tem-mean(tem, na.rm=TRUE)
        tem = cosine_taper(tem, taper = taper)
        tem =  c(tem,rep(0,krow-Ns))
        if(length(tem)<krow)
          {
            DMAT[,i] = rep(NA, length=krow)
          }
        else
          {
            DMAT[,i] = tem
          }
      }

    DFFT = mvfft(DMAT)

   DSPEC = Mod(DFFT)

    x = (ibeg+Ns/2)*dt

    freqs = df*c(0:((numfreqs/2)-1),(-numfreqs/2):(-1)  )

    y = (1:(numfreqs/2))*2*df

    RET = list(z=t(DSPEC[1:(numfreqs/2), ]), y=y, x=x, original_signal = sig, tt = seq(length(sig)) * dt)


    invisible(RET)
}
 

ftspec_image <- function(sig, dt, ft, time_span = NULL, freq_span = NULL, amp_span = NULL, taper = 0.05, scaling = "none", grid=TRUE, colorbar=TRUE, backcol=c(0, 0, 0), colormap=NULL, pretty=TRUE, ...)
{
	#Plots a Fourier spectrogram
	#INPUTS
        #    SIG is the signal to analyze
        #    DT is the sample rate (must be constant)
	#    FT is the Fourier transform input parameters, adopted from Jonathan Lees' code in RSEIS
	#        FT$NFFT is the fft length
	#        FT$NS is the number of samples in a window
        #        FT$NOV is the number of samples to overlap
        #    TIME_SPAN is the time span to plot, NULL plots everything
        #    FREQ_SPAN is the frequency span to plot (<=max frequency in spectrogram), NULL plots everything up to the Nyquist frequency
	#    AMP_SPAN is the amplitude range to plot.  NULL plots everything.
        #    TAPER is the cosine taper factor (amount of the signal to apply the taper to, must be < 0.5)
        #    SCALING determines whether to apply a natural log (ln), logarithmic (log), or square root (sqrt) scaling to the amplitude data
        #    GRID is a boolean asking whether to display grid lines
        #    COLORBAR is a boolean asking whether to plot an amplitude colorbar
        #    BACKCOL is a 3 element vector of RGB values for the background of the spectrogram, based on a 0 to 255 scale: [red, green, blue]
        #    COLORMAP is an R palette object determining how the spectrogram colors should look
        #    PRETTY is a boolean asking whether to adjust axis labels so that they're pretty (TRUE) or give the exactly specified time and frequency intervals (FALSE)
        #    OPTIONAL PARAMETERS
        #       TRACE_FORMAT is the format of the trace minima and maxima in sprintf format
        #       IMG_X_FORMAT is the format of the X axis labels of the image in sprintf format
        #       IMG_Y_FORMAT is the format of the Y axis labels of the image in sprintf format
        #       COLORBAR_FORMAT is the format of the colorbar labels in sprintf format   
        #       CEX.LAB is the font size of the image axis labels
        #       CEX.COLORBAR is the font size of the colorbar
        #       CEX.TRACE is the font size of the trace axis labels
        #       IMG_X_LAB is the X - axis label of the image, it defaults to "time"
        #       IMG_Y_LAB is the Y - axis label of the image, it defaults to "frequency" 
        #OUTPUTS
        #    IMG is the spectrogram	
	opts = list(...)

        if(!"img_x_lab" %in% names(opts))
        {
            opts$img_x_lab = "time"
        }

        if(!"img_y_lab" %in% names(opts))
        {
            opts$img_y_lab = "frequency"
        }

        if(is.null(time_span))
        {
                time_span=c(dt, length(sig) * dt)
        }

        if(time_span[2] > length(sig) * dt)
        {
                time_span[2]= length(sig) * dt
                warning("The requested spectrogram is longer than the actual signal.")
        }

        if(is.null(freq_span))
        {
                freq_span=c(0, 1/(dt * 2))
        }
        if(freq_span[2] > 1 / (dt * 2))
        {
                freq_span[2] = 1 / (dt * 2)
                warning("Requested maximum frequency is higher than the Nyquist frequency.")
        }

	ev=evolutive_fft(sig, dt, ft, freq_span, taper) #Calculate the Fourier spectrogram
        ev$tt = seq(length(sig)) * dt

        if(is.null(amp_span))
        {    
             amp_span = c(min(ev$z), max(ev$z))
        }

        img = list()
        img$z = ev$z[ev$x >= time_span[1] & ev$x <= time_span[2], ev$y >= freq_span[1] & ev$y <= freq_span[2]]
        img$x = ev$x[ev$x >= time_span[1] & ev$x <= time_span[2]]
        img$y = ev$y[ev$y >= freq_span[1] & ev$y <= freq_span[2]]

        img$z[img$z<amp_span[1]] = NA
        img$z[img$z>amp_span[2]] = amp_span[2]
        img$z[img$z == 0] = NA

        if(scaling == "ln") #Scale by natural log
        {
            img$z = log(img$z)
        }

        if(scaling == "log") #Log 10 scale
        {
            img$z = log10(img$z)
        }

        if(scaling == "sqrt") #Take the square root
        {
            img$z = sqrt(img$z)
        }

        trace = list()
        trace$sig = ev$original_signal[ev$tt >= time_span[1] & ev$tt <= time_span[2]]
        trace$tt = ev$tt[ev$tt >= time_span[1] & ev$tt <= time_span[2]]

        hht_package_plotter(img, trace, opts$img_x_lab, opts$img_y_lab, window = ft$ns / length(sig), colormap = colormap, backcol = backcol, pretty = pretty, grid = grid, colorbar = colorbar, opts = opts)
 
        invisible(img)

}		

hh_render <- function(hres, dt, dfreq, freq_span = NULL, time_span = NULL, scaling = "none", verbose = TRUE)
{
	#Renders a spectrogram of EMD or Ensemble EMD (EEMD) results.
	#INPUTS
	#	HRES is a matrix of data generated by EEMD_COMPILE or the output of HHTRANSFORM
	#	it represents a set on all time/frequency/amplitude points from the given EEMD run
        #       DT is the time resolution of the spectrogram.  Currently, if there is a hres$dt field, DT must be greater than or equal to hres$dt.
        #       this prevents subsample resolution.
        #       DFREQ is the frequency resolution of the spectrogram
        #       FREQ_SPAN is the frequency range to calculate the spectrum over c(MIN, MAX).  NULL means capture the full frequency spectrum of the signal.
        #       TIME_SPAN is the portion of the signal to include.  NULL means the whole signal.
        #       SCALING determines whether to plot frequency as log 10 ("log") or linear ("none")
        #       VERBOSE prints out status messages (i.e. IMF 1 COMPLETE!)
	#OUTPUTS
	#	HGRAM is a spectrogram matrix ready to be plotted by HHGRAM_IMAGE
	#Danny Bowman
	#UNC Chapel Hill

        hgram = hres

        if(scaling == "log")
        {
            hres$hinstfreq = log10(hres$hinstfreq)
        }
        else if (scaling != "none")
        {
             warning("Did not recognize scaling request \"", scaling, ".\" Reverting to linear frequency (scaling = \"none\").")
        }
       
        #Deal with logarithms of 0
        hres$hamp[hres$hinstfreq == -Inf] = 0
        hres$hinstfreq[hres$hinstfreq == -Inf] = 0

        if(is.null(freq_span))
        {
             freq_span = c(min(hres$hinstfreq), max(hres$hinstfreq))
        }
	
        if(!"trials" %in% names(hres))
 	{
		hres$trials=1
	}
        
        if("dt" %in% names(hres))
        {
             if(hres$dt > dt) #We don't want to have to interpolate between samples
             {
                 warning(paste("The time resolution", sprintf("%.2e", dt), "is lower than the sample rate", sprintf("%.2e", hres$dt), "of the time series.  This may introduce time gaps in the spectrogram."))
             }
             if("tt" %in% names(hres))
             {
                 warning("Input data has both DT (sample rate) and TT (sample times) components.  Component TT will be used to calculate the spectrogram")
                 hgram$tt = hres$tt
             }
             else
             {
                 hgram$tt = seq(length(hres$original_signal)) * hres$dt
             }
        }

        if(is.null(time_span))
        {
           time_span = c(min(hgram$tt), max(hgram$tt))
        }

        if(!(("tt" %in% names(hres)) | ("dt" %in% names(hres))))
        {
            warning("Neither DT (sample rate) nor TT (sample times) were specified in the input data.  Assuming DT is 1...")
            hgram$tt = seq(length(hres$original_signal))
        } 

        if(time_span[2]>max(hgram$tt))
        {
                time_span[2]=max(hgram$tt)
                warning("Requested time window is longer than the actual signal.")
        }
 
        hres$hinstfreq = array(hres$hinstfreq[which(hgram$tt >= time_span[1] & hgram$tt <= time_span[2]),], dim = c(length(hgram$tt), hres$nimf, hres$trials))
        hres$hamp = array(hres$hamp[hgram$tt >= time_span[1] & hgram$tt <= time_span[2],], dim = c(length(hgram$tt), hres$nimf, hres$trials))
        hres$original_signal = hres$original_signal[hgram$tt >= time_span[1] & hgram$tt <= time_span[2]]
        hgram$tt = hgram$tt[hgram$tt >= time_span[1] & hgram$tt <= time_span[2]]
       
        grid = list()
        grid$x = hgram$tt
        grid$y = seq(from = freq_span[1], to = freq_span[2] + dfreq, by = dfreq)
        hgram$z=array(rep(0,(length(grid$x) * length(grid$y) * hres$nimf)),dim=c(length(grid$x),length(grid$y), hres$nimf))
        hgram$cluster=hgram$z #Shows how many times a given grid node has data.
        for(i in seq(hres$nimf))
        {
            x = array(c(rep(hgram$tt,hres$trials), hres$hinstfreq[,i,]), dim = c(length(hgram$tt)*hres$trials, 2))
            imf_img = as.image(hres$hamp[,i,], grid = grid, x = x)
            hgram$z[,,i] = imf_img$z
            hgram$cluster[,,i] = imf_img$weights
            if(verbose)
            {
                print(paste("IMF", i, "COMPLETE!"))
            }
        }
       
        hgram$hinstfreq = hres$hinstfreq
        hgram$hamp = hres$hamp 
        hgram$z[is.na(hgram$z)] = 0
        hgram$cluster[is.na(hgram$cluster)] = 0
        hgram$x = imf_img$x 
        hgram$y = imf_img$y
	hgram$dfreq=dfreq
	hgram$dt=hres$dt
        hgram$scaling = scaling
	invisible(hgram) #Return the spectrogram structure.
}

hhgram_image <- function(hgram,time_span = NULL,freq_span = NULL, amp_span = NULL, clusterspec = FALSE, cluster_span=NULL, imf_list = NULL, imf_sum = FALSE, scaling = "none", grid=TRUE, colorbar=TRUE, backcol=c(0, 0, 0), colormap=NULL, pretty=FALSE, ...)
{
	#Plots a spectrogram of the EEMD processed signal as an image.	
	#INPUTS
	#	HGRAM is the subsetted spectrogram  from HH_RENDER.
	#		HGRAM$X is time
	#		HGRAM$Y is frequency
	#		HGRAM$Z is amplitude normalized to trials
	#		HGRAM$CLUSTER is a matrix containing integer values corresponding to the number of times a signal was recorded in a given spectrogram cell during EEMD
	#		The more often the signal is recorded, the more likely it is that the signal is real and not noise
	#		HGRAM$TRIALS is the number of times EEMD was run to generate signal
	#		HGRAM$ORIGINAL_SIGNAL is the original seismogram (without added noise)
        #               HGRAM$TT is the sample times
	#	TIME_SPAN is the time span to plot, [0,-1] plots everything
	#	FREQ_SPAN is the frequency span to plot (<=max frequency in spectrogram), [0,-1] plots everything
	#	AMP_SPAN is the amplitude span to plot, everything below is set to black, everything above is set to max color, [0, -1] scales to range in signal
        #	CLUSTERSPEC tells the code to plot the signal amplitude (FALSE) or the number of times data occupies a given pixel (TRUE).
	#	CLUSTER_SPAN plots only the parts of the signal that have a certain number of data points per pixel [AT LEAST, AT MOST] this only applies to EEMD with multiple trials.
        #       IMF_LIST is a list of IMFs to plot on the spectrogram.  If NULL, plot all IMFs.
        #       IMF_SUM can be set to show the sum of IMFs shown in the spectrogram plotted as a red line against the original trace
        #       SCALING determines whether to apply a logarithmic (log), or square root (sqrt) scaling to the amplitude data, default is "none"
	#	GRID is a boolean asking whether to display grid lines
	#	COLORBAR is a boolean asking whether to plot an amplitude colorbar
        #       BACKCOL is a 3 element vector of RGB values for the background of the spectrogram, based on a 0 to 255 scale: [red, green, blue]
        #       COLORMAP is an R palette object determining how the spectrogram colors should look
        #       PRETTY is a boolean asking whether to adjust axis labels so that they're pretty (TRUE) or give the exactly specified time and frequency intervals (FALSE)
	#OPTIONAL PARAMETERS
        #       TRACE_FORMAT is the format of the trace minima and maxima in sprintf format
        #       IMG_X_FORMAT is the format of the X axis labels of the image in sprintf format
        #       IMG_Y_FORMAT is the format of the Y axis labels of the image in sprintf format
        #       COLORBAR_FORMAT is the format of the colorbar labels in sprintf format   
        #       CEX.LAB is the font size of the image axis labels
        #       CEX.COLORBAR is the font size of the colorbar
        #       CEX.TRACE is the font size of the trace axis labels
        #       IMG_X_LAB is the X - axis label of the image, it defaults to "time"
        #       IMG_Y_LAB is the Y - axis label of the image, it defaults to "frequency"
        #OUTPUTS
        #     IMG is the spectrogram returned as an image

        opts = list(...)
 
        if(!"img_x_lab" %in% names(opts))
        {
            opts$img_x_lab = "time"
        }
        
        if(!"img_y_lab" %in% names(opts))
        {
            opts$img_y_lab = "frequency"
        }
  
        #Subset by IMFs
        if(is.null(imf_list))
        {
            imf_list = seq(hgram$nimf)
        }
        else
        {
            if(max(imf_list) > hgram$nimf)
            {
                warning("Requested more IMFs than are present in the actual EMD results!")
                imf_list = imf_list[imf_list < hgram$nimf]
            }
        }   
   
        if(!is.null(cluster_span))
        {
            img$z[cluster >= cluster_span[1] & cluster <= cluster_span[2]] = NA
        } 
 
	if(is.null(time_span))
	{
		time_span=c(min(hgram$tt), max(hgram$tt))
	}
	
	if(time_span[2]>max(hgram$tt))
	{
		time_span[2]=max(hgram$tt)
		warning("Requested time window is longer than the actual signal.")
	}
	
	if(is.null(freq_span))
	{
		freq_span=c(min(hgram$y), max(hgram$y))
	}
	if(freq_span[2]>max(hgram$hinstfreq))
	{
		freq_span[2]=max(hgram$y)
		warning("Requested frequency window is higher than maximum frequency in the spectrogram.")
	}

        if(imf_sum)
        {
             imf_sum = rowSums(hgram$averaged_imfs[hgram$x >= time_span[1] & hgram$x <= time_span[2], imf_list])
        }
        else
        {
            imf_sum = NULL
        }

        img = list()
        img$x = hgram$x[hgram$x >= time_span[1] & hgram$x <= time_span[2]]
        img$y = hgram$y[hgram$y >= freq_span[1] & hgram$y <= freq_span[2]]
        cluster = apply(hgram$cluster[hgram$x >= time_span[1] & hgram$x <= time_span[2], hgram$y >= freq_span[1] & hgram$y <= freq_span[2],imf_list], c(1, 2), sum)

        #Determine if we are plotting clustering or amplitudes

        if(clusterspec)
        {
            img$z = cluster
        }
        else
        {
            img$z = apply(hgram$z[hgram$x >= time_span[1] & hgram$x <= time_span[2], hgram$y >= freq_span[1] & hgram$y <= freq_span[2],imf_list], c(1, 2), sum)
        }

        if(!is.null(cluster_span))
        {
            img$z[cluster <= cluster_span[1] |  cluster >= cluster_span[2]] = NA
        }
       

        if(is.null(amp_span))
        {
             amp_span = c(min(img$z), max(img$z))
        }
 
        img$z[img$z<amp_span[1]] = NA
        img$z[img$z>amp_span[2]] = amp_span[2]
        img$z[img$z == 0] = NA

        if(scaling == "log") #Log 10 scale
        {
            img$z = log10(img$z) 
        }
   
        if(scaling == "sqrt") #Take the square root
        {
            img$z = sqrt(img$z)
        }
       
        trace = list()
        trace$sig = hgram$original_signal[hgram$tt >= time_span[1] & hgram$tt <= time_span[2]]
        trace$tt = hgram$tt[hgram$tt >= time_span[1] & hgram$tt <= time_span[2]]

        hht_package_plotter(img, trace, opts$img_x_lab, opts$img_y_lab, imf_sum = imf_sum, colormap = colormap, backcol = backcol, pretty = pretty, grid = grid, colorbar = colorbar, opts = opts)
    
        invisible(img)
}

hht_package_plotter <- function(img, trace, img_x_lab, img_y_lab, imf_sum = NULL, window = NULL, colormap = NULL, backcol = c(0, 0, 0), pretty = FALSE, grid = TRUE, colorbar = TRUE, opts = list())
{
    #Plots images and time series for Hilbert spectra, Fourier spectra, and cluster analysis.
    #This function is internal to the package and users should not be calling it.
    #
    #INPUTS
    #    IMG is the image portion of the figure
    #        IMG$X is the columns
    #        IMG$Y is the rows
    #        IMG$Z is the image data
    #    TRACE is the time series to plot at the top of the figure
    #        TRACE$SIG is the time series
    #        TRACE$TT is the time of each sample
    #    IMG_X_LAB is the label of the X axis of the image
    #    IMG_Y_LAB is the label of the Y axis of the image
    #    IMF_SUM is a red line on the time series plot showing the sum of the plotted IMFs, if available
    #        IMF_SUM$SIG is the summed IMFS
    #        IMF_SUM$TT is the time of each sample.  We assume all IMFS have equivalent timing.
    #    WINDOW is the length of the Fourier window, if applicable
    #    COLORMAP is the colormap to use for the image
    #    BACKCOL is the background color of the image
    #    PRETTY allows for nice axis labels
    #    GRID draws a grid on the image
    #    COLORBAR puts a colorbar corresponding to the range of values on the image
    #
    #    OPTS    OTHER POSSIBLE OPTIONS
    #        OPTS$TRACE_FORMAT is the format of the trace minima and maxima in sprintf format
    #        OPTS$IMG_X_FORMAT is the format of the X axis labels of the image in sprintf format
    #        OPTS$IMG_Y_FORMAT is the format of the Y axis labels of the image in sprintf format
    #        OPTS$COLORBAR_FORMAT is the format of the colorbar labels in sprintf format   
    #        OPTS$CEX.LAB is the font size of the image axis labels
    #        OPTS$CEX.COLORBAR is the font size of the colorbar
    #        OPTS$CEX.TRACE is the font size of the trace axis labels
    #        OPTS$TRACE_COL is the color of the trace
    #        OPTS$IMF_SUM_COL is the color of the IMF sums (if shown)

    #Configure parameters
    
    if(!"trace_format" %in% names(opts))
    {
        opts$trace_format = "%.1e"
    }
 
    if(!"img_x_format" %in% names(opts))
    {
        opts$img_x_format = "%.2f"
    }
   
    if(!"img_y_format" %in% names(opts))
    {
        opts$img_y_format = "%.2f"  
    }
  
    if(!"colorbar_format" %in% names(opts))
    {
         opts$colorbar_format = "%.1e"
    }
 
    if(!"cex.main" %in% names(opts))
    {
        opts$cex.main = 1
    }
    
    if(!"cex.trace" %in% names(opts))
    {
        opts$cex.trace = opts$cex.main * 0.75
    }

    if(!"cex.colorbar" %in% names(opts))
    {
        opts$cex.colorbar = opts$cex.main * 0.75
    }

    if(!"cex.lab" %in% names(opts))
    {
        opts$cex.lab = opts$cex.main
    }

    if(!"imf_sum_col" %in% names(opts))
    {
        opts$imf_sum_col = "red"
    }
   
    if(!"trace_col" %in% names(opts))
    {
        opts$trace_col = "black"
    }

    if(pretty)
    {   #Get nice divisions
        pretty_x = pretty(img$x, n=10)
        pretty_y = pretty(img$y, n=5) 
        pretty_x = pretty_x[pretty_x <= max(img$x) & pretty_x >= min(img$x)]
        pretty_y = pretty_y[pretty_y <= max(img$y) & pretty_y >= min(img$y)]
        img$z = img$z[img$x <= max(pretty_x) & img$x >= min(pretty_x), img$y <= max(pretty_y) & img$y >= min(pretty_y)]
        img$x = img$x[img$x <= max(pretty_x) & img$x >= min(pretty_x)]
        img$y = img$y[img$y <= max(pretty_y) & img$y >= min(pretty_y)]
        img_x_labels=sprintf(opts$img_x_format, pretty_x)
        img_y_labels=sprintf(opts$img_y_format, pretty_y)
        cat("Adjusting Time and Frequency limits to nice looking numbers (the \"pretty\" option is currently set to TRUE)\n")
    }    
    else 
    {    
       img_x_labels=sprintf(opts$img_x_format, seq(min(img$x), max(img$x), length.out = 10))
       img_y_labels=sprintf(opts$img_y_format, seq(min(img$y), max(img$y), length.out=5))
    }    

    if(is.null(colormap))
    {
        colormap=rainbow(500,start=0,end=5/6)
    }

    colorbins = length(colormap)
    
    plot(c(-0.15,1),c(-0.15,1),type="n",axes=FALSE,xlab="", ylab="") # Set up main plot window
 
    #Plot TRACE

    sig = trace$sig - mean(trace$sig)
    trace_y=0.75
    trace_x=0
    trace_yspan=0.10
    trace_xspan=0.9
    trace_at=seq(trace_y,trace_y+trace_yspan,length.out=2)
    trace_labels=c(min(trace$sig), max(trace$sig))
    trace_scale=trace_yspan/(max(sig)-min(sig))
    tt_scale=trace_xspan/(max(trace$tt) - min(trace$tt))
    axis(4,pos=trace_x+trace_xspan,at=trace_at, labels=c("",""), cex.axis=opts$cex.trace)
    lines((trace$tt - min(trace$tt)) * tt_scale + trace_x, trace_y + (sig - min(sig)) * trace_scale, col = opts$trace_col)
    if(!is.null(imf_sum))
    {
         lines(((trace$tt - min(trace$tt))*tt_scale+trace_x), (trace_y + (imf_sum - min(sig)) * trace_scale), col = opts$imf_sum_col)
    }
    rect(trace_x, trace_y, trace_x+trace_xspan, trace_y+trace_yspan)

    #Plot IMAGE
    
    image_y=0
    image_x=0
    image_yspan=0.75
    image_xspan=0.9
    image_xvec=seq(image_x, image_x+image_xspan, length.out=length(img$x))
    image_yvec=seq(image_y, image_y+image_yspan, length.out=length(img$y))
    img_x_at=seq(image_x,image_x+image_xspan,length.out=length(img_x_labels))
    img_y_at=seq(image_y,image_y+image_yspan, length.out=length(img_y_labels))
    rect(image_x,image_y,image_x+image_xspan,image_y+image_yspan,col=rgb(red=backcol[1], green=backcol[2], blue=backcol[3], maxColorValue=255))
    image(image_xvec,image_yvec, img$z, col=colormap,add=TRUE)
    axis(2, pos=image_x, at=img_y_at,labels=img_y_labels, cex.axis=opts$cex.lab)
    axis(1,pos=image_y, at=img_x_at,labels=img_x_labels, cex.axis=opts$cex.lab)

    #Plot Fourier window, if applicable
    
    if(!is.null(window))
    {
        rwidth = trace_xspan * window 
        rect(trace_x + trace_xspan - rwidth, trace_y + trace_yspan, trace_x + trace_xspan, trace_y + trace_yspan + 0.01, col = "black")
    }


    #Plot GRID
    if(grid)
    {
        line_color=rgb(red=100, green=100, blue=100, maxColorValue=255)
        line_type=3
        for(k in 2:(length(img_x_at)-1))
        {
            lines(c(img_x_at[k], img_x_at[k]), c(image_y, trace_y+trace_yspan), col=line_color, lty=line_type)
        }

        for(k in 2:(length(img_y_at)-1))
        {
            lines(c(image_x, image_x+image_xspan), c(img_y_at[k], img_y_at[k]), col=line_color, lty=line_type)
        }
    }

    #Plot COLORBAR

    if(colorbar)
    {
        color_x=image_x+image_xspan+0.015
        color_xspan=0.025
        color_y=image_y+image_yspan-0.20
        color_yspan=0.10
        color_xvec=c(color_x,color_x+color_xspan)
        color_yvec=seq(color_y, color_y+color_yspan, length.out=colorbins)
        color_at=seq(color_y,color_y+color_yspan,length.out=2)
        colorbar_matrix=array(seq_len(colorbins),dim=c(1, colorbins))
        image(color_xvec, color_yvec, colorbar_matrix, col=colormap, axes=FALSE, add=TRUE)
    }


    #Plot TEXT


    text(trace_x + trace_xspan + 0.03, trace_y, srt=90, sprintf(opts$trace_format,trace_labels[1]), cex=opts$cex.trace)
    text(trace_x + trace_xspan + 0.03, trace_y+trace_yspan, srt=90, sprintf(opts$trace_format, trace_labels[2]), cex=opts$cex.trace)
    text(image_x-0.095, image_y+image_yspan/2, srt=90, img_y_lab, cex=opts$cex.lab)
    text(image_x+image_xspan/2, image_y-0.1, img_x_lab, cex=opts$cex.lab)
    if("main" %in% names(opts))
    {
        text(trace_x+trace_xspan/2, trace_y+trace_yspan+0.05,opts$main, cex=opts$cex.main)
    }
    if(colorbar)
    {
        text(color_x+0.015, color_y-0.0125, sprintf(opts$colorbar_format, min(img$z[!is.na(img$z)])), cex=opts$cex.colorbar)
        text(color_x+0.015, color_y+color_yspan+0.0125, sprintf(opts$colorbar_format, max(img$z[!is.na(img$z)])), cex=opts$cex.colorbar)
    }
 
}

hh_spectrum <- function(hres, dt, dfreq, freq_span = NULL, time_span = NULL, scaling = "none", verbose = TRUE)
{
    #Calculate the Hilbert spectrogram of a signal contained in HRES (returned by HHTRANSFORM or EEMD_COMPILE)
    #INPUTS
    #       HRES is a matrix of data generated by EEMD_COMPILE or the output of HHTRANSFORM
    #       it represents a set on all time/frequency/amplitude points from the given EEMD run
    #       DT is the time resolution of the spectrogram.  Currently, if there is a hres$dt field, DT must be greater than or equal to hres$dt.
    #       this prevents subsample resolution.
    #       DFREQ is the frequency resolution of the spectrogram
    #       FREQ_SPAN is the frequency range to calculate the spectrum over c(MIN, MAX).  NULL means capture the full frequency spectrum of the signal.
    #       TIME_SPAN is the time span to calculate the spectrum over c(MIN, MAX).  NULL means use the entire signal
    #       SCALING determines whether to calculate frequency as log 10 ("log") or linear ("none")
    #       VERBOSE prints out status messages (i.e. IMF 1 COMPLETE!)
    #OUTPUTS
    #    HSPEC is the Hilbert spectrum of the signal, separated by IMF.

   hgram = hh_render(hres, dt, dfreq, freq_span = NULL, time_span = NULL, scaling = scaling, verbose = TRUE)

   amps = array(rep(0, dim(hgram$z)[2] * dim(hgram$z)[3]), dim = dim(hgram$z)[2:3])

   for(i in seq(hres$nimf))
   {
       amps[, i] = apply(hgram$z[, , i], 2, sum) 
   }

  invisible(list(amplitude = amps, frequency = hgram$y, original_signal = hgram$original_signal, dt = hgram$dt))
} 

hhspec_plot <- function(hspec, freq_span = NULL, scaling = "none", imf_list = NULL, show_total = TRUE, show_fourier = FALSE, show_imfs = FALSE, legend = TRUE, ...)
{
    #Plot the Hilbert spectrum, optionally as individual IMFs, optionally with the scaled Fourier spectrum for comparison
    #INPUTS
    #    HSPEC is the Hilbert spectrogram returned by HHSPECTRUM
    #    FREQ_SPAN is the frequencies to plot, NULL means plot everything
    #    SCALING whether to take the base 10 logarithm of amplitude ("log") or square root of amplitude ("sqrt")  or do nothing ("none")
    #    IMF_LIST means only include these IMFS, NULL includes all of them
    #    SHOW_TOTAL means show the sum of the IMF Hilbert spectra
    #    SHOW_IMFS means plot individual IMFs
    #    SHOW_FOURIER determines whether you want a Fourier spectrum for comparison (TRUE) or not (FALSE)
    #    IMF_COLS is a vector of length IMF_LIST with colors to plot the individual IMFs.  Defaults to a colormap
    #    LEGEND asks whether to plot a legend.  Additional options will place the legend where you want it.
    #ADDITIONAL OPTIONS
    #    XLAB is the X axis label
    #    YLAB is the Y axis label
    #    LEGEND_LOCATION determines where to put the legend.
    #    TOTAL_COL is the color of the ensemble Hilbert spectrum
    #    TOTAL_LWD is the line weight of the ensemble Hilbert spectrogram
    #    LOTAL_LTY is the line type of the ensemble Hilbert spectrogram
    #    IMF_COLS sets the color of each IMF - a vector with length IMF_LIST    
    #    IMF_LWD is the line weight for the IMFs as a vector with length IMF_LIST
    #    IMF_LTY is the line type for the IMFs as a vector with length IMF_LIST
    #    FOURIER_COL is the color of the Fourier spectrum line
    #    FOURIER_LTY is the line type of the Fourier spectrum line
    #    FOURIER_LWD is the line weight of the Fourier spectrum line
    #    SCALE_FOURIER scales the Fourier spectrum line to the Hilbert spectrum line if TRUE.  Defaults to FALSE.

    if(!(show_total | show_imfs | show_fourier))
    {
        error("Nothing to plot!  Set at least one of SHOW_TOTAL, SHOW_IMFS, or SHOW_FOURIER to TRUE.")
    }

    opts = list(...)

    if(!(scaling == "log" | scaling == "sqrt" | scaling == "none"))
    {
        warning(paste("Did not recognize requested scaling: \"", scaling, "\".  Options are \"log\" (base 10 logarithm), \"sqrt\" (square root), or \"none\""))
        scaling = "none"
    }
    
    if(is.null(freq_span))
    {
        freq_span = c(0, max(hspec$frequency))
    }
   
    hspec$amplitude = hspec$amplitude[hspec$frequency >= freq_span[1] & hspec$frequency<= freq_span[2],]
    hspec$frequency = hspec$frequency[hspec$frequency >= freq_span[1] & hspec$frequency<= freq_span[2]]

    if(!"legend_location" %in% names(opts) & legend)
    {
        opts$legend_location = "topright"
    }


    if(!"total_col" %in% names(opts))
    {
        opts$total_col = "red"
    }

    if(!"total_lwd" %in% names(opts))
    {
        opts$total_lwd = 1
    }
    
    if(!"total_lty" %in% names(opts))
    {
        opts$total_lty = 1
    }

    if(!"xlab" %in% names(opts))
    {
        opts$xlab = "frequency"
    }

    if(!"ylab" %in% names(opts))
    {
        if(scaling != "none")
        {
            opts$ylab = paste(scaling, "amplitude")
        }
        else
        {
             opts$ylab = "amplitude"
        }
    }
    
    if(is.null(imf_list))
    {
        imf_list = seq(dim(hspec$amplitude)[2])
    }

    if(!"imf_cols" %in% names(opts))
    {
        if(show_total)
        {
            opts$imf_cols = rainbow(length(imf_list), start = 1/6, end = 5/6)
        }
        else
        {
            opts$imf_cols = rainbow(length(imf_list), start = 0, end = 5/6)
        }
    }
   
    if(!"imf_lwd" %in% names(opts))
    {
        opts$imf_lwd = rep(1, length(imf_list))
    }

    if(!"imf_lty" %in% names(opts))
    {
        opts$imf_lty = rep(1, length(imf_list))
    }

    if(!"fourier_col" %in% names(opts))
    {
        opts$fourier_col = "black"
    }

    if(!"fourier_lty" %in% names(opts))
    {
        opts$fourier_lty = 1
    }
   
    if(!"fourier_lwd" %in% names(opts))
    {
        opts$fourier_lwd = 1
    }

    if(!"scale_fourier" %in% names(opts))
    {
        opts$scale_fourier = FALSE
    }
    
    pmin = Inf
    pmax = -Inf

    if(show_imfs)
    {
        imf_amp = hspec$amplitude[, imf_list]
        pmin = min(imf_amp[imf_amp>0])
        pmax = max(imf_amp)
    }

    if(show_total)
    {
        if(length(imf_list)>1)
        {
            total_amp = apply(hspec$amplitude[,imf_list], 1, sum)
        }
        else
        {
            total_amp = hspec$amplitude[,imf_list]
        }
        if(max(total_amp) > pmax)
        {
            pmax = max(total_amp[total_amp > 0])
        }
        if(min(total_amp) < pmin)
        {
            pmin = min(total_amp[total_amp > 0])
        }
    }

     if(show_fourier)
     {
        fourier_freqs = seq(0, 1/(hspec$dt * 2), length.out = length(hspec$original_signal)-1)
        fspec = Mod(fft(hspec$original_signal - mean(hspec$original_signal)))[1:length(hspec$original_signal)/2][fourier_freqs >= freq_span[1] & fourier_freqs <= freq_span[2]]
        if(opts$scale_fourier)
        {
            fspec = fspec * pmax/max(fspec)
        }
        if(max(fspec) > pmax)
        {
            pmax = max(fspec)
        }
        if(min(fspec[fspec > 0]) < pmin)
        {
            pmin = min(fspec[fspec > 0])
        }
    }
    
    if(scaling == "log")
    {
        pmax = log10(pmax)
        pmin = log10(pmin)
    }

    if(scaling == "sqrt")
    {
        pmax = sqrt(pmax)
        pmin = sqrt(pmin)
    }
    
    plot(c(min(hspec$frequency), max(hspec$frequency)), c(pmin, pmax), type = "n", xlab = opts$xlab, ylab = opts$ylab)

    if(show_imfs)
    {
       for(k in seq(length(imf_list)))
       {

           amp = imf_amp[,k]
           if(scaling == "log")
           {
              amp = log10(amp)
           }

           if(scaling == "sqrt")
           {
               amp = sqrt(amp)
           } 
           lines(hspec$frequency[amp > -Inf], amp[amp > -Inf], col = opts$imf_cols[k], lwd = opts$imf_lwd[k], lty = opts$imf_lty[k])
       }
    }
   
    if(show_total) 
    {
        if(scaling == "log")
        {
            total_amp = log10(total_amp)
        }

        if(scaling == "sqrt")
        {
            total_amp = sqrt(total_amp)
        }

        lines(hspec$frequency, total_amp, lwd = opts$total_lwd, lty = opts$total_lty, col = opts$total_col)
    }

    if(show_fourier)
    {

        if(scaling == "log")
        {
            fspec = log10(fspec)
        }

        if(scaling == "sqrt")
        {
            fspec = sqrt(fspec)
        }

        lines(fourier_freqs[fourier_freqs >= freq_span[1] & fourier_freqs <= freq_span[2]], fspec, 
            lty = opts$fourier_lty, lwd = opts$fourier_lwd, col = opts$fourier_col)
    }

    if(legend)
    {
        legend_labs = c()
        legend_cols = c()
        legend_lty = c()
        legend_lwd = c()
        if(show_total)
        {
            legend_labs = c(legend_labs, "Total Hilbert")
            legend_cols = c(legend_cols, opts$total_col)
            legend_lty = c(legend_lty, opts$total_lty)
            legend_lwd = c(legend_lwd, opts$total_lwd) 
        }
        if(show_imfs) 
        {
            legend_labs = c(legend_labs, paste(rep("IMF", length(imf_list)), imf_list))
            legend_cols = c(legend_cols, opts$imf_cols)
            legend_lty = c(legend_lty, opts$imf_lty)
            legend_lwd = c(legend_lwd, opts$imf_lwd)
        }

        if(show_fourier)
        {
            legend_labs = c(legend_labs, "Fourier")
            legend_cols = c(legend_cols, opts$fourier_col)
            legend_lty = c(legend_lty, opts$fourier_lty[1])
            legend_lwd = c(legend_lwd, opts$fourier_lwd[1])
        }
        legend(opts$legend_location, legend = legend_labs, lty = legend_lty, lwd = legend_lwd, col = legend_cols)
     }
}

central_difference <- function(sig, tt)
{
    #Calculate derivatives using the central difference method for a (possibly irregular) time series
    #INPUTS
    #    SIG is the time series
    #    TT is the sample times, increasing only
    #OUTPUTS
    #    DSIG is the derivative of SIG

    lx=length(sig)
    f2=sig[3:lx]
    t2 = tt[3:lx]
    f1 = sig[1:(lx-2)]
    t1 = tt[1:(lx - 2)]

    #central difference
    df0=(f2-f1)/(t2 - t1)

    #forward difference for the first point
    df1=(sig[2]-sig[1])/(tt[2] -tt[1])

    #backward difference for the last point
    df2=(sig[lx]-sig[lx-1])/(tt[lx] - tt[lx - 1])

    dsig=c(df1,df0,df2)
    invisible(dsig)
}

hilbert_transform <- function(sig)
{
   #Return the Hilbert transform of a signal.
   #Code modified from the EMD package by Donghoh Kim and Hee-Seok Oh (http://dasan.sejong.ac.kr/~dhkim/software_emd.html)
   #INPUTS
   #    SIG - the signal to be transformed
   #OUTPUTS
   #    ASIG - the analytic_signal
   ndata = length(sig)
   h = rep(0, ndata)

   if(ndata %% 2 == 0)
   {
       h[c(1, ndata/2+1)] = 1 
       h[2:(ndata/2)] = 2 
   }
   else
   {
       h[1] = 1
       h[2:((ndata + 1)/2)] = 2 
   }

   asig = fft(h * fft(sig), inverse = TRUE)/ndata
   invisible(asig)
} 

instantaneous_frequency <- function(asig, tt)
{
    #Calculate instantaneous frequency via method outlined in Equation 6 of
    #Dasios, A.; Astin, T. R. & McCann, C. Compressional-wave Q estimation from full-waveform sonic data 
    #Geophysical Prospecting, 2001, 49, 353-373
    #INPUTS
    #    ASIG is the analytic signal 
    #    TT is the sample times
    #OUTPUTS
    #    INSTFREQ is the instantaneous frequency 

    dsig = central_difference(Re(asig), tt)
    dhsig = central_difference(Im(asig), tt)

    instfreq = (1/(2*pi))*(Re(asig)*dhsig - Im(asig)*dsig)/((Re(asig)^2) + (Im(asig)^2))
    
    invisible(instfreq)
}
    
hilbert_envelope <- function(asig)
{
    #Calculate the envelope (instantaneous amplitude) of a signal.
    #INPUTS
    #    ASIG is the analytic signal
    #OUTPUTS
    #    ENVELOPE is the positive envelope of the signal

    envelope = abs(asig)
    invisible(envelope)
}
        
plot_imfs <-function(sig, time_span, imf_list, original_signal, residue, fit_line=FALSE, lwd=1, cex=1, ...)
{
    #Better IMF plotter
    #This function plots IMFs on the same plot so they can be checked for mode mixing or other problems.
    #It plots shifted traces in a single window
    #INPUTS
    #    SIG is the signal data structure returned by EEMD or EMD analysis
    #    Note that SIG$AVERAGED_IMFS will be plotted instead of SIG$IMF, likewise SIG$AVERAGED_RESIDUE takes precidence
    #    over SIG$RESIDUE, if both exist.
    #        SIG$IMF is a N by M array where N is the length of the signal and M is the number of IMFs
    #        SIG$ORIGINAL_SIGNAL is the original signal before EEMD
    #        SIG$RESIDUE is the residual after EMD
    #        SIG$DT is the sample rate
    #    TIME_SPAN is a 2 element vector giving the time range to plot
    #    IMF_LIST is the IMFs to plot
    #    ORIGINAL_SIGNAL is a boolean asking if you are going to plot the original signal also (defaults to be on top)
    #    RESIDUE is a boolean asking if you are going to plot the residual (defaults to be on bottom)
    #	 FIT_LINE is a boolean asking if you want to plot a line showing the sum of IMFs on top of the original signal (to check how the selected IMFs reconstruct the original signal)
    #	 LWT is the line weight (for plotting figures)
    #    CEX is the size of text (for plotting figures)
    #    ... other parameters to pass to main plotting function
   
 
    if(time_span[2]<0)
    {
        time_span[2]=length(sig$original_signal)*sig$dt
    }
    
    if(time_span[1]==0)
    {
        time_span[1]=sig$dt
    }
    
    if("averaged_imfs" %in% names(sig))
    {
        sig$imf=sig$averaged_imfs
    }

    if("averaged_residue" %in% names(sig))
    {
        sig$residue=sig$averaged_residue
    }


    time_ind=ceiling(time_span[1]/sig$dt):floor(time_span[2]/sig$dt)
    tt=time_ind*sig$dt
    
    plot(c(0, 1), c(0, 1), type="n", axes=FALSE, xlab="Time (s)", ylab="", cex.lab=cex, ...)
    
    yax_labs=c()
    snum=length(imf_list)+residue+original_signal
    sp=1/snum # Spacing of subplots

    if(original_signal)
    {
         snum=snum+1
         os=sig$original_signal[time_ind]-mean(sig$original_signal[time_ind])
         scale=max(abs(os)) 
    }
    else
    {
        scale=max(abs(sig$imf))
    }
    
    if(residue)
    {
        snum=snum+1
        res=sig$residue[time_ind]-mean(sig$residue[time_ind])
	res=res*(sp/(2*scale))
        yax_labs=append(yax_labs,"Residue")
    }

    
    trace_pos=sp/2 #Where the trace starts on the plot
    imfs=sig$imf*(sp/(scale*2))
    ts=(tt-min(tt))*(1/(time_span[2]-time_span[1]))

    if(residue)
    {
        lines(ts, res+trace_pos, lwd=lwd)
        trace_pos=trace_pos+sp
    } 
    
    for(k in rev(imf_list))
    {
       lines(ts, imfs[time_ind,k]+trace_pos, lwd=lwd)
       trace_pos=trace_pos+sp
       yax_labs=append(yax_labs, paste("IMF",k))
    }
    if(original_signal)
    {
        lines(ts, os*(sp/(2*scale))+trace_pos, lwd=lwd)
        yax_labs=append(yax_labs,"Signal")
        if(fit_line)
        {
            if(length(imf_list)>1)
            {
                fline=rowSums(imfs[time_ind,imf_list])
            }
            else
            {
                fline=imfs[time_ind,imf_list]
            }
            if(residue)
            {
                fline=fline+res
            }
            lines(ts, fline+trace_pos, lwd=lwd, col="red")
        }
    }
    xax_labs=pretty(seq(min(tt)-sig$dt, max(tt), length.out=11))
    axis(1, pos=0, at=seq(0,1, length.out=length(xax_labs)), labels=xax_labs, cex.axis=cex)
    axis(2, pos=0, at=seq(sp/2, 1, by=sp), labels=yax_labs, cex.axis=cex)
    segments(c(0,0,1, 0), c(0, 1, 1, 0), c(0,1, 1, 1), c(1,1, 0, 0), lwd=lwd) 
}

sig2imf <- function(sig, tt, emd_config = list(tol = 5, max_sift = 200, stop_rule = "type5", boundary = "wave", sm = "none", smlevels = c(1), spar = NULL, max_imf = 100, interm = NULL))

{
    #Extract IMFs
    #This function is intended to take data, recover IMFs, and save them
    #It calls and runs code developed by Donghoh Kim and Hee-Seok Oh as part of the "EMD" package available on CRAN.
    #Refer to their documentation for details on the EMD algorithm as implemented here.
    #INPUTS
    #	SIG is the time series
    #   TT is the sample times 
    #	EMD_CONFIG controls how the EMD algorithm operates. This documentation adapted from the "emd" function documentation in the EMD package, current as of June 2013
    #         EMD_CONFIG$MAX_SIFT stop sifting after this many times
    #         EMD_CONFIG$STOP_RULE as quoted from the EMD package:  "stopping rule of sifting. The type1 stopping rule indicates that absolute values 
    #         of envelope mean must be less than the user-specified tolerance level in the sense
    #         that the local average of upper and lower envelope is zero. The stopping rules
    #         type2, type3, type4 and type5 are the stopping rules given by equation (5.5)
    #         of Huang et al. (1998), equation (11a), equation (11b) and S stoppage of Huang
    #         and Wu (2008), respectively."
    #         EMD_CONFIG$BOUNDARY - how the beginning and end of the signal are handled
    #         EMD_CONFIG$SM - Specifies how the signal envelope is constructed, see Kim et al, 2012.
    #         EMD_CONFIG$SMLEVELS - Specifies what level of the IMF is obtained by smoothing other than interpolation - not sure what this means
    #         EMD_CONFIG$SPAR - User-defined smoothing parameter for spline, kernal, or local polynomial smoothign
    #         EMD_CONFIG$MAX_IMF - How many IMFs are allowed, IMFs above this number will not be recorded
    #         EMD_CONFIG$INTERM - specifies vector of periods to be excluded from IMFs to cope with mode mixing.  I do not use this; instead I use the EEMD method.
    #OUTPUT is a list containing the original signal, IMFs, and information on EMD parameters.
    #Danny Bowman
    #UNC Chapel Hill
    
    #NOTES ON STOP RULES
    #These are rules to determine when sifting should stop and the result be saved as an IMF
    #The rules are implemented in the dcb_extractimf, which is taken almost verbatim from the extractimf function in the EMD package.
    #Here is how they work
    #TYPE1
    #type 1 finds the mean of the envelopes as defined by the spline function.  See Fig. 3b in Huang et al 1998:
    #The envelope is the dotted lines, the mean is the solid line.  In an IMF, the authors of the EMD function
    #think that this mean should be close to 0 (hence, below the defined tolerance).  This will force the IMF to be
    #symmetric, but I am concerned it may err towards frequency modulation at the expense of amplitude.
    #This criterion is not mentioned in literature I have read and they only cite Huang et al 1998 where
    #it is certainly not mentioned at all.  Instead they have derived this criterion from the definition of an IMF:
    #that at every point the mean defined by the local maximum and minimum envelope be zero.
    #TYPE2
    #Type 2 is the method used by Huang et al 1998.  It finds the standard deviation between two consecutive siftings;
    #when the SD falls below a certain tolerance value the sifting stops.  Huang et al recommend a tolerance value of 0.2 or 0.3
    #TYPE3
    #Type 3 is my addition to the EMD package.  It implements the S stoppage criterion per Huang et al 2008 (further described in
    #Huang et al 2003).  The S criterion means the sifting process stops when there has been no change in the number of extrema
    #and zero crossings for S iterations.  S is arbitrary but Huang's tests indicate around 3-8 works well but this should be tested
    #per Huang et al 2003 for each type of data otherwise it is somewhat arbitrary.

    emd_result=emd(sig, tt, max.sift=emd_config$max_sift, stoprule=emd_config$stop_rule, tol=emd_config$tol, 
        boundary=emd_config$boundary,sm=emd_config$sm,spar=emd_config$spar, 
        check=FALSE, plot.imf=FALSE,max.imf=emd_config$max_imf)
    emd_result$original_signal=sig
    emd_result$tt=tt
    for(pname in names(emd_config))
    {
        emd_result[pname]=emd_config[pname]
    }
  
    emd_result$hinstfreq = array(rep(0, length(emd_result$original_signal) * emd_result$nimf), dim = c(length(emd_result$original_signal), emd_result$nimf))
    emd_result$hamp = emd_result$hinstfreq
    
    for(i in seq(emd_result$nimf))
    {
        imf = emd_result$imf[,i]
        aimf = hilbert_transform(imf)
        emd_result$hinstfreq[, i] = instantaneous_frequency(aimf, tt)
        emd_result$hamp[, i] = hilbert_envelope(aimf)
    }
    invisible(emd_result)
}

precision_tester <- function(tt = seq(0, 10, by = 0.01), a = 1, b = 1, c = 1, omega_1 = 2 * pi, omega_2 = 4 * pi, phi_1 = 0, phi_2 = pi/6, plot_signal = TRUE, plot_instfreq = TRUE, plot_error = TRUE, ...)
{
    #This function computes the instantaeous frequency of a signal of the form
    # a sin(omega_1 t + phi_1) + b sin(omega_2 + phi_2) + c
    #where a, b, c, omega_1, omega_2, phi_1, and phi_2 are real numbers.
    #The instantaneous frequency is calculated in two ways:

    #1.  An exact analytic expression calulated for signals of this form using Equation 6 in 
    #Dasios, A.; Astin, T. R. & McCann, C. Compressional-wave Q estimation from full-waveform sonic data 
    #Geophysical Prospecting, 2001, 49, 353-373.
    #An exact expression can be derived through the liberal use of algebra and trigonometric identities, see http://www.unc.edu/~haksaeng/hht/analytic_instantaneous_freq.pdf

    #2.  Using the numeric method presented in this R package (i.e. the functions hilbert_transform and instantaneous_frequency)
    
    #The PRECISION_TESTER function allows a comparison of these two methods - my hope is it may identify where the numeric method falls short.

    #INPUTS
    #    TT - sample times
    #    A - Amplitude coefficient of first sinusoid
    #    B - Amplitude coefficient of second sinusoid
    #    C - Constant shift
    #    OMEGA_1 - frequency of first sinusoid, in radians
    #    OMEGA_2 - frequency of second sinusoid, in radians
    #    PHI_1 - phase shift of first sinusoid, in radians
    #    PHI_2 - phase shift of second sinusoid, in radians
    #    PLOT_SIGNAL - If TRUE, show the sinusoid defined by the above parameters
    #    PLOT_INSTFREQ - If TRUE, plot the analytic and numeric instantaneous frequencies against each other
    #    PLOT_ERROR - If TRUE, plot the error between the analytic and numeric instantaneous frequencies
    #    ... passes plot parameters to plotter
    #OUTPUTS
    #    INSTFREQ is the instantaneous frequency and the time series
    #        INSTFREQ$SIG is the  time series defined by the input parameters
    #        INSTFREQ$ANALYTIC is the analytically calculated frequency
    #        INSTFREQ$NUMERIC is the frequency calculated via this package's numeric algorithm

    A = sin(omega_1 * tt + phi_1)
    B = sin(omega_2 * tt + phi_2)
    C = cos(omega_1 * tt + phi_1)
    D = cos(omega_2 * tt + phi_2)
   
    #Time series
    sig = a * A + b * B  + c

    #Instantaneous frequency derived analytically

    num = omega_1 * a^2 + omega_2 * b^2 + a * b * (omega_1 + omega_2) * (A * B + C * D) + c * (omega_1 * a * A + omega_2 * b * B)
    denom = a^2 + 2 * a * b * (A * B + C * D) + 2 * c * (a * A + b * B) + b^2 + c^2

    analytic_instfreq = num / (denom * 2 * pi)

    #Instantaneous frequency derived numerically

    asig = hilbert_transform(sig)

    numeric_instfreq = instantaneous_frequency(asig, tt)

    #PLOTTING 

    if(plot_signal)
    {
        dev.new()
        plot(tt, sig, type = "l", xlab = "Time", ylab = "Amplitude", main = "Time series", ...)
    }

    if(plot_instfreq)
    {
        dev.new()
        ylow = min(c(min(analytic_instfreq), min(numeric_instfreq)))
        yhigh = max(c(max(analytic_instfreq), max(numeric_instfreq)))
        plot(tt, analytic_instfreq, type = "l", col = "red", ylim = c(ylow, yhigh), xlab = "Time", ylab = "Frequency", main = "Analytically and numerically derived values for instantaneous frequency", ...)
        points(tt, numeric_instfreq, ...)
        legend("topright", lty = c(1, NA), pch = c(NA, 1), legend = c("Analytic", "Numeric"), col = c("black", "red"))
    }

    if(plot_error)
    {
        dev.new()
        plot(tt, analytic_instfreq - numeric_instfreq, type = "l", xlab = "Time", ylab = "Frequency Error", main = "Numerically derived instantaneous frequency subtracted from analytically derived instantaneous frequency", ...)
    }

    instfreq = list(sig = sig, analytic = analytic_instfreq, numeric = numeric_instfreq)
    invisible(instfreq)
}




