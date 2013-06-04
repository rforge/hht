## THIS COLLECTION OF FUNCTIONS IMPLEMENTS EMD VIA THE "EMD" PACKAGE
##AND ALSO RUNS THE ENSEMBLE EMD METHOD

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

sig2imf <- function(sig, tt, tol = 5, max_sift = 200, stop_rule = "type5", boundary = "wave", sm = "none", smlevels = c(1), spar = NULL, max_imf = 100, interm = NULL)

{
    #Extract IMFs
    #This function is intended to take data, recover IMFs, and save them
    #It calls and runs code developed by Donghoh Kim and Hee-Seok Oh as part of the "EMD" package available on CRAN.
    #Refer to their documentation for details on the EMD algorithm as implemented here.
    #INPUTS
    #	SIG is the time series
    #   TT is the sample times 
    #   MAX_SIFT stop sifting after this many times
    #   STOP_RULE as quoted from the EMD package:  "stopping rule of sifting. The type1 stopping rule indicates that absolute values 
    #   of envelope mean must be less than the user-specified tolerance level in the sense
    #   that the local average of upper and lower envelope is zero. The stopping rules
    #   type2, type3, type4 and type5 are the stopping rules given by equation (5.5)
    #   of Huang et al. (1998), equation (11a), equation (11b) and S stoppage of Huang
    #   and Wu (2008), respectively."
    #   BOUNDARY - how the beginning and end of the signal are handled
    #   SM - Specifies how the signal envelope is constructed, see Kim et al, 2012.
    #   SMLEVELS - Specifies what level of the IMF is obtained by smoothing other than interpolation - not sure what this means
    #   SPAR - User-defined smoothing parameter for spline, kernal, or local polynomial smoothign
    #   MAX_IMF - How many IMFs are allowed, IMFs above this number will not be recorded
    #   INTERM - specifies vector of periods to be excluded from IMFs to cope with mode mixing.  I do not use this; instead I use the EEMD method.
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

    emd_result=emd(sig, tt, max.sift=max_sift, stoprule=stop_rule, tol=tol, 
        boundary=boundary,sm=sm,spar=spar, 
        check=FALSE, plot.imf=FALSE,max.imf=max_imf)
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
