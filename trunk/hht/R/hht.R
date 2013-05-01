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
		hht=hhtransform(emd_result)
		emd_result$hinstfreq=hht$hinstfreq
		emd_result$hamp=hht$hamp
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
	averaged_imfs=averaged_imfs/counter
	averaged_noise=averaged_noise/counter
	averaged_residue=averaged_residue/counter
	
        if(counter<trials)
        {
                warning("Number of trials requested was greater than the number of trials found in the trials directory")
        }

	eemd_result=c()
        eemd_result$nimf=nimf
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

    if(krow - Ns < 0) #I was having problems with this earlier.  The number of rows became a negative number.  Not sure how you want to handle this.
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

    RET = list(sig=sig, dt=dt, numfreqs=numfreqs, wpars=list(Nfft= Nfft,  Ns=Ns, Nov=Nov, fl=freq_span[1], fh=freq_span[2]), DSPEC=DSPEC, freqs=y, tims=x)


    invisible(RET)
}
 

ftspec_image <- function(sig, dt, ft, time_span, freq_span, amp_span, amp_units=NULL, amp_unit_conversion=NULL, taper = 0.05, grid=TRUE, colorbar=TRUE, backcol=c(0, 0, 0), colormap=NULL, pretty=TRUE, cex=1, trace_format = "%1.1e", main="")
{
	#Plots a Fourier spectrogram
	#INPUTS
        #SIG is the signal to analyze
        #DT is the sample rate
	#FT is the Fourier transform input parameters, adopted from Jonathan Lees' code in RSEIS
	#	FT$NFFT is the fft length
	#	FT$NS is the number of samples in a window
	#	FT$NOV is the number of samples to overlap
        #TIME_SPAN is the time span to plot, [0,-1] plots everything
        #FREQ_SPAN is the frequency span to plot (<=max frequency in spectrogram), [0,-1] plots everything up to the Nyquist frequency
	#AMP_SPAN is the amplitude range to plot.  [0, -1] plots everything.
        #AMP_UNITS is the units of amplitude, used for axes labels
        #AMP_UNIT_CONVERSION tells how to convert units in input data to units to display
        #GRID is a boolean asking whether to display grid lines
        #COLORBAR is a boolean asking whether to plot an amplitude colorbar
        #BACKCOL is a 3 element vector of RGB values for the background of the spectrogram, based on a 0 to 255 scale: [red, green, blue]
        #COLORMAP is an R palette object determining how the spectrogram colors should look
        #PRETTY is a boolean asking whether to adjust axis labels so that they're pretty (TRUE) or give the exactly specified time and frequency intervals (FALSE)
        #CEX is the size of the text on the figure
        #TRACE_FORMAT is the number format of the trace labels
        #MAIN is the title of the plot.
	
	
	options(digits.secs=3)

        if(time_span[2]<0)
        {
                time_span[2]=length(sig)*dt
        }

        if(time_span[2]>length(sig)*dt)
        {
                time_span[2]=length(sig)*dt
                warning("Requested time window is longer than the actual signal.")
        }

        if(freq_span[2]<0)
        {
                freq_span[2]=(1/dt)/2
        }
        if(freq_span[2]>(1/dt)/2)
        {
                freq_span[2]=(1/dt)/2
                warning("Requested frequency window is higher than the Nyquist frequency.")
        }

        if(is.null(amp_unit_conversion))
        {
                amp_unit_conversion=1
        }
        
        if(is.null(colormap))
        {
            colormap=rainbow(500,start=0,end=5/6)
        }

        colorbins=length(colormap)

        if(pretty)
        {
            time_labels=pretty(seq(0,time_span[2]-time_span[1],length.out=10)+time_span[1], n=10)
            freq_labels=pretty(seq(freq_span[1],freq_span[2],length.out=5), n=5)
            time_span=c(min(time_labels), max(time_labels))
            freq_span=c(min(freq_labels), max(freq_labels))
        }
        else
        {
           time_labels=format(seq(0,time_span[2]-time_span[1],length.out=10)+time_span[1],digits=2)
           freq_labels=format(seq(freq_span[1],freq_span[2],length.out=5), digits=2)
        }

        tind=seq_len((time_span[2]-time_span[1])/dt)+(time_span[1])/dt

        par(mai=c(0.1, 0.1, 0.1, 0.1))
        plot(c(-0.15,1),c(-0.15,1),type="n",axes=FALSE,xlab="", ylab="") # Set up opts$main plot window
        par(mai=c(0.1, 0.1, 0.1, 0.1))

        sig=sig[tind]*amp_unit_conversion


	ev=evolutive_fft(sig, dt, ft, freq_span)

        ampcex=0.5*cex

        #Plot TRACE

        tt=tind*dt-time_span[1]
        trace_y=0.75
        trace_x=0
        trace_yspan=0.10
        trace_xspan=0.9
        trace_at=seq(trace_y,trace_y+trace_yspan,length.out=2)
        trace_labels=c(min(sig), max(sig))
        trace_scale=trace_yspan/(max(sig)-min(sig))
        tt_scale=trace_xspan/max(tt)
        axis(4,pos=trace_x+trace_xspan,at=trace_at, labels=c("",""), cex.axis=ampcex)
        lines((tt*tt_scale+trace_x), (trace_y+(sig-min(sig))*trace_scale))
        rect(trace_x, trace_y, trace_x+trace_xspan, trace_y+trace_yspan)

        #Plot IMAGE
	
        f_ind=(ev$freqs>=freq_span[1] & ev$freqs <= freq_span[2])
        freqs=ev$freqs[f_ind]
        image_z=t(ev$DSPEC[1:(ev$numfreqs/2),])
        image_z=image_z[,f_ind]
        if(amp_span[2]<0)
        {

             amp_span[1]=min(image_z)
             amp_span[2]=max(image_z)
        }
        image_z[image_z<amp_span[1]]=NA
        image_z[image_z>amp_span[2]]=amp_span[2]
        image_y=0
        image_x=0
        image_yspan=0.75
        image_xspan=0.9
        image_xvec=seq(image_x, image_x+image_xspan, length.out=length(ev$tims))
        image_yvec=seq(image_y, image_y+image_yspan, length.out=length(freqs))
        time_at=seq(image_x,image_x+image_xspan,length.out=length(time_labels))
        freq_at=seq(image_y,image_y+image_yspan, length.out=length(freq_labels))
        rect(image_x,image_y,image_x+image_xspan,image_y+image_yspan,col=rgb(red=backcol[1], green=backcol[2], blue=backcol[3], maxColorValue=255))
        image(image_xvec,image_yvec, image_z, col=colormap,add=TRUE)
        axis(2, pos=image_x, at=freq_at,labels=freq_labels, cex.axis=cex)
        axis(1, pos=image_y, at=time_at,labels=time_labels, cex.axis=cex)
        
        #Plot GRID
        
        if(grid)
        {
                line_color=rgb(red=100, green=100, blue=100, maxColorValue=255)
                line_type=3
                for(k in 2:(length(time_at)-1))
                {
                        lines(c(time_at[k], time_at[k]), c(image_y, trace_y+trace_yspan), col=line_color, lty=line_type)
                }

                for(k in 2:(length(freq_at)-1))
                {
                        lines(c(image_x, image_x+image_xspan), c(freq_at[k], freq_at[k]), col=line_color, lty=line_type)
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
                color_labels=format(seq(amp_span[1],amp_span[2],length.out=2),digits=2)
                colorbar=array(seq_len(colorbins),dim=c(1, colorbins))
                image(color_xvec, color_yvec, colorbar, col=colormap, axes=FALSE, add=TRUE)
                text(color_x+0.015, color_y-0.0125, format(amp_span[1], digits=2), cex=ampcex)
                text(color_x+0.015, color_y+color_yspan+0.0125, format(amp_span[2], digits=2), cex=ampcex)        
        }

        #Plot TEXT
        text(trace_x + trace_xspan + 0.03, trace_y, srt=90, sprintf(opts$trace_format,trace_labels[1]), cex=ampcex)
        text(trace_x + trace_xspan + 0.03, trace_y+trace_yspan, srt=90, sprintf(opts$trace_format, trace_labels[2]), cex=ampcex)
	text(image_x-0.095, image_y+image_yspan/2, srt=90, "Frequency", cex=cex)
        text(image_x+image_xspan/2, image_y-0.1, "Time (s)", cex=cex)
        text(trace_x+trace_xspan/2, trace_y+trace_yspan+0.05,opts$main, cex=cex)
        
        if(!is.null(amp_units))
        {
                text(trace_x+trace_xspan/2, trace_y+trace_yspan+0.025,paste("Trace and Spectrogram Amplitudes in",amp_units), cex=ampcex)
        }
 
       #Plot WINDOW
       
       rwidth=trace_xspan*(ev$wpars$Ns/length(sig))		
       rect(trace_x+trace_xspan-rwidth, trace_y+trace_yspan, trace_x+trace_xspan, trace_y+trace_yspan+0.01, col="black")

}		

hh_render <- function(hres, dt, dfreq, freq_span = NULL, verbose = TRUE)
{
	#Renders a spectrogram of EMD or Ensemble EMD (EEMD) results.
	#INPUTS
	#	HRES is a matrix of data generated by EEMD_COMPILE or the output of HHTRANSFORM
	#	it represents a set on all time/frequency/amplitude points from the given EEMD run
        #       DT is the time resolution of the spectrogram.  Currently, if there is a hres$dt field, DT must be greater than or equal to hres$dt.
        #       this prevents subsample resolution.
        #       DFREQ is the frequency resolution of the spectrogram
        #       FREQ_SPAN is the frequency range to calculate the spectrum over c(MIN, MAX).  NULL means capture the full frequency spectrum of the signal.
        #       VERBOSE prints out status messages (i.e. IMF 1 COMPLETE!)
	#OUTPUTS
	#	HSPEC is a spectrogram matrix ready to be plotted by HHSPEC_IMAGE
	#Danny Bowman
	#UNC Chapel Hill
	#August 2012

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
             }
        }

	#If the spectra is generated from EMD instead of EEMD, make it resemble EEMD but with only 1 trial.
	if(length(dim(hres$hamp))==2)
	{
		hres$hamp=array(hres$hamp, dim=c(dim(hres$hamp)[1], dim(hres$hamp)[2], 1))
		hres$hinstfreq=array(hres$hinstfreq, dim=c(dim(hres$hamp)[1], dim(hres$hamp)[2], 1))
	}
      
        if(!(("tt" %in% names(hres)) | ("dt" %in% names(hres))))
        {
            warning("Neither DT (sample rate) nor TT (sample times) were specified in the input data.  Assuming DT is 1...")
            hres$dt = 1
        } 

        if(!"tt" %in% names(hres))
        { 
            tt = seq(0, dim(hres$hamp)[1] - 1) * hres$dt #Time of signal
        }
        else #Possible irregular sampling
        {
            tt = hres$tt
        }
        
        hspec = hres
        grid = list()
        grid$x = tt
        grid$y = seq(from = freq_span[1], to = freq_span[2] + dfreq, by = dfreq)
        hspec$z=array(rep(0,(length(grid$x) * length(grid$y) * hres$nimf)),dim=c(length(grid$x),length(grid$y), hres$nimf))
        hspec$cluster=hspec$z #Shows how many times a given grid node has data.
        for(i in seq(hres$nimf))
        {
            x = array(c(rep(tt,hres$trials), hres$hinstfreq[,i,]), dim = c(length(tt)*hres$trials, 2))
            imf_img = as.image(hres$hamp[,i,], grid = grid, x = x)
            hspec$z[,,i] = imf_img$z
            hspec$cluster[,,i] = imf_img$weights
            if(verbose)
            {
                print(paste("IMF", i, "COMPLETE!"))
            }
        }
        
        hspec$z[is.na(hspec$z)] = 0
        hspec$cluster[is.na(hspec$cluster)] = 0
        hspec$x = imf_img$x 
        hspec$y = imf_img$y
	
	hspec$original_signal=hres$original_signal
	hspec$dfreq=dfreq
	hspec$dt=hres$dt
        hspec$tt = hres$tt
	invisible(hspec) #Return the spectrogram structure.
}

hhspec_image <- function(hspec,time_span,freq_span, amp_span, cluster_span=NULL, imf_list = NULL, imf_trace = FALSE, scaling = "none", grid=TRUE, colorbar=TRUE, backcol=c(0, 0, 0), colormap=NULL, pretty=FALSE, ...)
{
	#Plots a spectrogram of the EEMD processed signal as an image.	
	#INPUTS
	#	HSPEC is the subsetted spectrogram  from HH_RENDER.
	#		HSPEC$X is time
	#		HSPEC$Y is frequency
	#		HSPEC$Z is amplitude normalized to trials
	#		HSPEC$CLUSTER is a matrix containing integer values corresponding to the number of times a signal was recorded in a given spectrogram cell during EEMD
	#		The more often the signal is recorded, the more likely it is that the signal is real and not noise
	#		HSPEC$TRIALS is the number of times EEMD was run to generate signal
	#		HSPEC$ORIGINAL_SIGNAL is the original seismogram (without added noise)
        #               HSPEC$TT is the sample times
	#	TIME_SPAN is the time span to plot, [0,-1] plots everything
	#	FREQ_SPAN is the frequency span to plot (<=max frequency in spectrogram), [0,-1] plots everything
	#	AMP_SPAN is the amplitude span to plot, everything below is set to black, everything above is set to max color, [0, -1] scales to range in signal
	#	CLUSTER_SPAN plots only the parts of the signal that have a certain number of data points per pixel [AT LEAST, AT MOST] this only applies to EEMD with multiple trials.
        #       IMF_LIST is a list of IMFs to plot on the spectrogram.  If NULL, plot all IMFs.
        #       IMF_TRACE can be set to show the sum of IMFs shown in the spectrogram plotted as a red line against the original trace
        #       SCALING determines whether to apply a natural log (ln), logarithmic (log), or square root (sqrt) scaling to the amplitude data
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

	#OUTPUTS
	#	HHIMAGE is the result of the requested spectrogram subsetting.


        opts = list(...)

        #Subset by IMFs
        if(is.null(imf_list))
        {
            imf_list = seq(hspec$nimf)
        }
        else
        {
            if(max(imf_list) > hspec$nimf)
            {
                warning("Requested more IMFs than are present in the actual EMD results!")
                imf_list = imf_list[imf_list < hspec$nimf]
            }
        }   
   
        img = list()
        img$z = apply(hspec$z[,,imf_list], c(1, 2), sum)
        img$x = hspec$x
        img$y = hspec$y
        cluster = apply(hspec$cluster[,,imf_list], c(1, 2), sum)
        
        if(!is.null(cluster_span))
        {
            img$z[cluster >= cluster_span[1] & cluster <= cluster_span[2]] = NA
        } 
 
	if(time_span[2]<0)
	{
		time_span[2]=max(hspec$tt)
	}
	
	if(time_span[2]>max(hspec$tt)
	{
		time_span[2]=max(hspec$tt)
		warning("Requested time window is longer than the actual signal.")
	}
	
	if(freq_span[2]<0)
	{
		freq_span[2]=max(hspec$hinstfreq)
	}
	if(freq_span[2]>hspec$max_freq)
	{
		freq_span[2]=hspec$max_freq
		warning("Requested frequency window is higher than maximum frequency in the spectrogram.")
	}

        imf_sum = rowSums(hspec$averaged_imfs[hspec$x >= time_span[1] & hspec$x <= time_span[2], imf_list])
        img = list()
        img$z = apply(hspec$z[hspec$x >= time_span[1] & hspec$x <= time_span[2], hspec$y >= freq_span[1] & hspec$y <= freq_span[2],imf_list], c(1, 2), sum)
        img$x = hspec$x[hspec$x >= time_span[1] & hspec$x <= time_span[2]]
        img$y = hspec$y[hspec$y >= freq_span[1] & hspec$y <= freq_span[2]]
        cluster = apply(hspec$cluster[hspec$x >= time_span[1] & hspec$x <= time_span[2], hspec$y >= freq_span[1] & hspec$y <= freq_span[2],imf_list], c(1, 2), sum)

        if(!is.null(cluster_span))
        {
            img$z[cluster <= cluster_span[1] |  cluster >= cluster_span[2]] = NA
        }

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
        
        hht_package_plotter(img, hspec$original_signal, img_x_lab, img_y_lab, imf_sum, colormap = colormap, backcol = backcol, pretty, grid, colorbar, ...)
}

hht_package_plotter <- function(img, trace, img_x_lab, img_y_lab, imf_sum = NULL, window = NULL, colormap = NULL, backcol = c(0, 0, 0), pretty = FALSE, grid = TRUE, colorbar = TRUE, opts)
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
    #    ... also passes lots of other possible options, such as the label of the axes, etc
    #
    #OPTS    OTHER POSSIBLE OPTIONS
    #    TRACE_FORMAT is the format of the trace minima and maxima in sprintf format
    #    IMG_X_FORMAT is the format of the X axis labels of the image in sprintf format
    #    IMG_Y_FORMAT is the format of the Y axis labels of the image in sprintf format
    #    COLORBAR_FORMAT is the format of the colorbar labels in sprintf format   
    #    CEX.LAB is the font size of the image axis labels
    #    CEX.COLORBAR is the font size of the colorbar
    #    CEX.TRACE is the font size of the trace axis labels

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

    if(pretty)
    {    
        img_x_labels=sprintf(opts$img_x_format, pretty(img$x, n=10))
        img_y_labels=sprintf(opts$img_y_format, pretty(img$y, n=5))
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

    trace_y=0.75
    trace_x=0
    trace_yspan=0.10
    trace_xspan=0.9
    trace_at=seq(trace_y,trace_y+trace_yspan,length.out=2)
    trace_labels=c(min(trace$sig), max(trace$sig))
    trace_scale=trace_yspan/(max(trace$sig)-min(trace$sig))
    tt_scale=trace_xspan/max(trace$tt)
    axis(4,pos=trace_x+trace_xspan,at=trace_at, labels=c("",""), cex.axis=opts$cex.trace)
    lines((trace$tt*tt_scale+trace_x), (trace_y+(trace$sig-min(trace$sig))*trace_scale))
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
    text(image_x-0.095, image_y+image_yspan/2, srt=90, img_x_lab, cex=opts$cex.lab)
    text(image_x+image_xspan/2, image_y-0.1, img_y_lab, cex=opts$cex.lab)
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


hhtransform <- function(imf_set)

{
    #Transform IMFs into instantaneous frequency and amplitude using the Hilbert transform
    #INPUTS
    #EMD_RESULT is the EMD decomposition of a signal returned by sig2imf.R
    #Danny Bowman
    #OUTPUTS
    #  HHT is the Hilbert Transform of EMD_RESULT
    #UNC Chapel Hill
   
    
   if("averaged_imfs" %in% names(imf_set))
    {
        imf_set$imf=imf_set$averaged_imfs
    }

    tt=seq(from=0, by=imf_set$dt, length=length(imf_set$original_signal))
    hht_result=imf_set
    hilbert=hilbertspec(imf_set$imf,tt=tt)
    hht_result$hamp=hilbert$amplitude
    hht_result$hinstfreq=hilbert$instantfreq
    invisible(hht_result)
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

sig2imf <- function(sig, dt, emd_config)

{
    #Extract IMFs
    #This function is intended to take data, recover IMFs, and save them
    #It calls and runs code developed by Donghoh Kim, Hee-Seok Oh as part of the "EMD" package available on CRAN.
    #I have modified their code and included some of it in this repository.
    #INPUTS
    #	SIG is the time series
    #   DT is the sample rate
    #	EMD_CONFIG controls how the EMD algorithm operations
    #         EMD_CONFIG$MAX_SIFT how many times the IMFs can be sifted
    #         EMD_CONFIG$MAX_IMF maximum number of IMFs that can be generated
    #         EMD_CONFIG$TOL Sifting stop criterion.
    #         EMD_CONFIG$STOP_RULE Make sure to read section on stop rules and make sure you understand what they imply!
    #         EMD_CONFIG$BOUNDARY How the start and stop of the time series are handled duing the spline process.
    #         EMD_CONFIG$SM Spline smoothing
    #         EMD_CONFIG$SPAR Smoothing parameter (only needed if sm is not none)
    #         EMD_CONFIG$WEIGHT Weight if "sm" is "spline"
    #         VERBOSE lets you know how many IMFs have been extracted.

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

    tt=seq_len(length(sig))*dt
    tt=tt[which(!is.na(sig))]
    sig=sig[which(!is.na(sig))]
    emd_result=emd(sig, tt, max.sift=emd_config$max_sift, stoprule=emd_config$stop_rule, tol=emd_config$tol, 
        boundary=emd_config$boundary,sm=emd_config$sm,spar=emd_config$spar,weight=emd_config$weight, 
        check=FALSE, plot.imf=FALSE,max.imf=emd_config$max_imf)
    emd_result$original_signal=sig
    emd_result$dt=dt
    for(pname in names(emd_config))
    {
        emd_result[pname]=emd_config[pname]
    }
    invisible(emd_result)
}
