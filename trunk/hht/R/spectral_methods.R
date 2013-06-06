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
        tem = spec.taper(tem, p = taper)
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

instantaneous_frequency <- function(asig, tt, method = "arctan", lag = 1)
{
    #Calculate instantaneous frequency via method outlined in Equation 6 of
    #Dasios, A.; Astin, T. R. & McCann, C. Compressional-wave Q estimation from full-waveform sonic data 
    #Geophysical Prospecting, 2001, 49, 353-373
    #INPUTS
    #    ASIG is the analytic signal 
    #    TT is the sample times
    #    METHOD is the way the differentiation is performed. "arctan" uses the arctangent of the real and imaginary parts of the Hilbert transform, taking the derivative of phase for frequency
    #    "chain" uses the analytical derivative of the arctangent function prior to performing the calculation 
    #    One must be cautious when using these two methods - they should be tested using precision_tester on signals with similar frequency and sample rate
    #    LAG determines the order of differentiation.  The default is 2 (central difference method), but this may be unstable at frequencies approaching the Nyquist frequency.
    #OUTPUTS
    #    INSTFREQ is the instantaneous frequency 

    if(!method %in% c("arctan", "chain"))
    {
        stop(paste("Did not recognize frequency calculation method:", method, "Please use either arctan or chain.", sep = " "))
    }

    if(method == "arctan") #This is the method used in the EMD package
    {
        phase = atan2(Im(asig), Re(asig))
        d = c(0, -diff(phase))
        p = 2 * pi  * ((d > pi) - (d < -pi))
        unphase = phase + cumsum (p)
        instfreq = abs(diff(unphase, lag) / diff(tt, lag))
        instfreq = abs(instfreq[-length(instfreq)] + instfreq[-1])/2 
        instfreq =  c(rep(instfreq[1], 1), instfreq, rep(instfreq[length(instfreq)], lag))
    }
    
    if(method == "chain") #Dasios et al 2001 "Compressional-wave Q estimation from full-waveform sonic data," Equation 6
    {
        dsig = diff(Re(asig), lag)/diff(tt, lag)
        dsig = c(rep(dsig[1], 1), dsig, rep(dsig[length(dsig)], lag - 1))
        dhsig = diff(Im(asig), lag)/diff(tt, lag)
        dhsig = c(rep(dhsig[1]), dhsig, rep(dhsig[length(dhsig)], lag - 1))
        instfreq = (Re(asig)*dhsig - Im(asig)*dsig)/((Re(asig)^2) + (Im(asig)^2))
    }
    
    invisible(instfreq/(2 * pi))
}
    

precision_tester <- function(tt = seq(0, 10, by = 0.01), method = "arctan", lag = 1, a = 1, b = 1, c = 1, omega_1 = 2 * pi, omega_2 = 4 * pi, phi_1 = 0, phi_2 = pi/6, plot_signal = TRUE, plot_instfreq = TRUE, plot_error = TRUE, ...)
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
    #    METHOD - either "arctan" or "chain", see function "instantaneous_frequency"
    #    LAG - Differentiation lag, see function "instantaneous frequency"
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

    numeric_instfreq = instantaneous_frequency(asig, tt, method = method, lag = lag)

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
        legend("topright", lty = c(1, NA), pch = c(NA, 1), legend = c("Analytic", "Numeric"), col = c("red", "black"))
    }

    if(plot_error)
    {
        dev.new()
        plot(tt, analytic_instfreq - numeric_instfreq, type = "l", xlab = "Time", ylab = "Frequency Error", main = "Numerically derived instantaneous frequency subtracted from analytically derived instantaneous frequency", ...)
    }

    instfreq = list(sig = sig, analytic = analytic_instfreq, numeric = numeric_instfreq)
    invisible(instfreq)
}

