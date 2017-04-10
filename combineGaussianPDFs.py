# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Scripts for assessing and combing measurements and upper/lower limits using Guassian PDFs
# E.g., obtaining a combined flux measurement with uncertainties from line fluxes upper limits
# measured for the same object from multiple spectra/observations.
#
# See 160828 notes for details
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import combineGaussianPDFs as cgp
import kbsutilities as kbs
import scipy.special
import scipy.interpolate as si
import scipy.signal
import numpy as np
import pdb
import sys
import matplotlib.pyplot as plt
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plotGauss1DRepresentation(values,errors,Nsigma=1,axislabel='line flux',Npoints=1000,
                              outputname='gaussvaluerep.pdf',skewparam=0,sigrange=5.0,
                              xpriorhat=False,legendlocation='upper left',usexlog=False,logxlow=1e-3,
                              magnifications=False,magerrors=False,fixxrange=False,verbose=True):
    """

    Plot measurements in 1D using Gaussians to represent the measureents and their uncertainties.
    Upper/lower lilmits should be given as fluxe measurements and uncertainty estimate as well.
    This gives a first of whether the independent measurements are in agreement with each other.
    The code will also return the best values from the combined profiles (product and sum) of
    these measurements and the corresponding standard deviation and sample 68% and 95% uncertainties.

    --- INPUT ---
    values          The actual measurements. For upper/lower limits a flux measurement should also be provided.
                    At the very least use the 1sigma uncertaintiy to indicate a 1sigma detection.
    errors          The estiamted uncertainteis on the measurements
    Nsigma          The number of sigma the errors represent
    axislabel       Label to put on x-axis of plot
    Npoints         Number of points to generate for each gaussian. Will be generate out to +/-sigrange sigma
    outputname      Name of output figure to stor plot to
    sewparam        To plot skewed Gaussians, provide list of skewness parameters
    sigrange        Range of sigmas to consider around gaussian means, i.e., mean+/-sigrange.
                    This is also the range considered to be 'within agreement' with each other, i.e., if
                    any of the measurements differ by more than sigrange it will drive the combined product to 0.
    xpriorhat       Use this keyword to put a top-hat prior on the gaussians by restricting the axxepted xvalues
                    to a range [xmingood,xmaxgood]. This is useful for, for instance restrcting fluxes to be
                    postive definit, by setting xpriorhat=[0,1e5]
    legendlocation  Providing the location of the legend within the plot. If False no legend will be plotted
    usexlog         To plot the x-axis in log set usexlog=True.
    logxlow         If usexlog=True use this keyword to provide the lower bound of the log-axes
    magnifications  Provide a list of lens magnifications if the measurements should be corrected for this.
    magerrors       The estimated uncertainty on the lensing magnifications.
    fixxrange       FOr instance for larger magnification errors the automatically set xrange for plottig
                    (set to include all gaussian representations) might not be ideal. To overwrite this for
                    plottin purposes provide the deisred x-range to this keyword. logxlow overwrite the lower bound.
    versose         Toggle verbosity

    --- OUTPUT ---
    sumvals         The best estimates and corresponding uncertainties estimated on the "PDF" from summing
                    the individual Gaussian representations. The output is a list:
                        sumvals  = [valmed,valmean,valstd,m2sig,m1sig,center,p1sig,p2sig]
                    contianing the following values:
                    valmed    Median of distribution sample
                    valmean   Mean of distribution sample
                    valstd    Standard deviation of distribution sample
                    m2sig     95% confindence interval (1sigma) lower bound around the central value
                    m1sig     68% confindence interval (1sigma) lower bound around the central value
                    center    Central value of distribution sample
                    p1sig     68% confindence interval (1sigma) upper bound around the central value
                    p2sig     95% confindence interval (2sigma) upper bound around the central value
    prodvals        The same as 'sumvals' but for the product "PDF" of the individual Gaussian representations.

    --- EXAMPLE OF USE ---
    import combineGaussianPDFs as cgp

    values = [0.7,2,2.2,4]
    errors = [1.2,1.3,1.25,1.5]
    sumvals, prodvals = cgp.plotGauss1DRepresentation(values,errors,axislabel='test',xpriorhat=[0,8])

    errors = [0.2,0.3,0.25,0.5]
    sumvals, prodvals = cgp.plotGauss1DRepresentation(values,errors,axislabel='test',sigrange=3)

    values  = [2,2]
    errors  = [1.2,1.2]
    mags    = [1,10]
    magerrs = [0,3]
    sumvals, prodvals = cgp.plotGauss1DRepresentation(values,errors,magnifications=mags,magerrors=magerrs,sigrange=3,xpriorhat=[0,8],outputname='gaussvaluerep_161115.pdf',usexlog=False,logxlow=1e-2,legendlocation='upper right')

    """
    values        = np.asarray(values) # make sure values is a numpy array
    errors        = np.asarray(errors) # make sure errors is a numpy array
    Nmeasurements = len(values)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Building dictionary with Guassian representations '
    gaussreps = {}

    # intitale xrange limits
    xmin = np.min(values)
    xmax = np.max(values)

    for vv, val in enumerate(values):
        xlow       = val-errors[vv]/Nsigma*sigrange
        xhigh      = val+errors[vv]/Nsigma*sigrange
        xvals_obs  = np.linspace(xlow,xhigh,Npoints)
        param      = val, errors[vv], skewparam
        if not magnifications:
            gaussvals  = cgp.gauss_skew(xvals_obs,*param)
            xvals      = xvals_obs
        else:
            if len(np.atleast_1d(magnifications)) == 1:
                magnifications = np.atleast_1d(magnifications).tolist()*Nmeasurements
            if len(np.atleast_1d(magerrors)) == 1:
                magerrors = np.atleast_1d(magerrors).tolist()*Nmeasurements

            xvals_true, gaussvals, gaussobs = cgp.correct4magnification(xvals_obs,val,errors[vv],
                                                                        magnifications[vv],magerrors[vv],
                                                                        Nsigma=sigrange,Npoints=Npoints)
            xvals = xvals_true

        if xpriorhat:
            gaussvals[xvals < xpriorhat[0]] = 0.0
            gaussvals[xvals > xpriorhat[1]] = 0.0

        gaussreps[str(vv)] = xvals, gaussvals
        xmin = np.min([xmin,np.min(xvals)])
        xmax = np.max([xmax,np.max(xvals)])

    xrange   = [xmin,xmax]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Generating combined profile of all gaussians '
    Npointscomb = Npoints
    xvals_comb  = np.linspace(xrange[0],xrange[1],Npointscomb)

    gaussarray  = np.zeros([Nmeasurements,Npointscomb])

    for vv, val in enumerate(values):
        xvals, gauss     = gaussreps[str(vv)]

        gauss_interp     = si.interp1d(xvals,gauss,bounds_error=False,fill_value=0.0)(xvals_comb)
        #gauss_interp     = kbs.interp(xvals,gauss,xvals_comb) # not set up to deal with xnew outside xold
        gaussarray[vv,:] = gauss_interp


    gauss_comb_sum   = np.sum(gaussarray,axis=0)
    gauss_comb_prod  = np.prod(gaussarray,axis=0)
    maxprod          = np.max(gauss_comb_prod)

    # gauss_comb_sum   = gauss_comb_sum/np.sum(gauss_comb_sum)
    # gauss_comb_prod  = gauss_comb_prod/np.sum(gauss_comb_prod)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Get measurement and uncertainty from sum- and product-combined PDFs'
    if verbose: print '   sampling "sum" distribution'
    Nsample       = 10000

    probabilities = gauss_comb_sum/np.sum(gauss_comb_sum)
    sum_sampled   = np.random.choice(xvals_comb, Nsample, p=probabilities)
    if verbose: print '   sampling "product" distribution'


    if verbose: print '   getting best values and uncertainties for output '
    valarr   = sum_sampled
    goodval  = valarr[(valarr != -99) & (valarr != 0)]
    valmed   = np.median(goodval)
    valmean  = np.mean(goodval)
    valstd   = np.std(goodval)
    Ngoodval = len(goodval)
    p2sig    = np.sort(goodval)[int(np.floor(0.975*Ngoodval))] # 95% conf upper bound
    p1sig    = np.sort(goodval)[int(np.floor(0.84*Ngoodval))]  # 68% conf upper bound
    center   = np.sort(goodval)[int(np.round(0.50*Ngoodval))]
    m1sig    = np.sort(goodval)[int(np.ceil(0.16*Ngoodval))]   # 68% conf lower bound
    m2sig    = np.sort(goodval)[int(np.ceil(0.025*Ngoodval))]  # 95% conf lower bound

    sumvals  = [valmed,valmean,valstd,m2sig,m1sig,center,p1sig,p2sig]


    if maxprod != 0.0:
        probabilities = gauss_comb_prod/np.sum(gauss_comb_prod)
        prod_sampled  = np.random.choice(xvals_comb, Nsample, p=probabilities)

        valarr   = prod_sampled
        goodval  = valarr[(valarr != -99) & (valarr != 0)]
        valmed   = np.median(goodval)
        valmean  = np.mean(goodval)
        valstd   = np.std(goodval)
        Ngoodval = len(goodval)
        p2sig    = np.sort(goodval)[int(np.floor(0.975*Ngoodval))] # 95% conf upper bound
        p1sig    = np.sort(goodval)[int(np.floor(0.84*Ngoodval))]  # 68% conf upper bound
        center   = np.sort(goodval)[int(np.round(0.50*Ngoodval))]
        m1sig    = np.sort(goodval)[int(np.ceil(0.16*Ngoodval))]   # 68% conf lower bound
        m2sig    = np.sort(goodval)[int(np.ceil(0.025*Ngoodval))]  # 95% conf lower bound

        prodvals = [valmed,valmean,valstd,m2sig,m1sig,center,p1sig,p2sig]
    else:
        prodvals = [0.0,0.0,-99.0,-99.0,-99.0,0.0,-99.0,-99.0]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up and generating plot '+outputname
    fig = plt.figure(figsize=(6, 3))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.98, bottom=0.13, top=0.96)
    Fsize    = 9
    lthick   = 2
    marksize = 4
    capsize  = 0
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)

    #-------------------------------------------------------
    if verbose: print '   plot sample histograms of sum and product distributions '
    Nbin = 200
    histvals,binvals,patches=plt.hist(sum_sampled,bins=Nbin,range=xrange,normed=True,
                                      facecolor='green',alpha=0.4,histtype='stepfilled',label='Sample from sum')

    if maxprod != 0:
        histvals,binvals,patches=plt.hist(prod_sampled,bins=Nbin,range=xrange,normed=True,
                                          facecolor='red',alpha=0.4,histtype='stepfilled',label='Sample from product')
        prodscale = np.max(gauss_comb_sum)/maxprod
        scalestr  = '\nscaled to max(sum)'
    else:
        prodscale = 1.0
        scalestr  = '\n(All 0s: No overlap for sigrange='+str(sigrange)+')'
    #-------------------------------------------------------
    if verbose: print '   plot sum and product distributions '
    plt.plot(xvals_comb,gauss_comb_sum,ls='-',color='green',lw=lthick*2,label='Sum of Gaussians')

    plt.plot(xvals_comb,gauss_comb_prod*prodscale,
             ls='-',color='red',lw=lthick*2,label='Product of Gaussians'+scalestr)

    #-------------------------------------------------------
    if verbose: print '   "best value" points above sum and product distributions'

    #[valmed,valmean,valstd,m2sig,m1sig,center,p1sig,p2sig]


    plt.errorbar(sumvals[0],np.max(gauss_comb_sum)*1.23,xerr=None,yerr=None,
                 fmt='*',lw=lthick, markersize=marksize,capsize=capsize,
                 markerfacecolor='green',ecolor='green',markeredgecolor='k',label='Median')
    plt.errorbar(sumvals[1],np.max(gauss_comb_sum)*1.18,xerr=sumvals[2],yerr=None,
                 fmt='D',lw=lthick, markersize=marksize,capsize=capsize,
                 markerfacecolor='green',ecolor='green',markeredgecolor='k',label='Mean $\pm$ STD')

    xerrlow  = np.asarray([sumvals[5]-sumvals[4]])
    xerrhigh = np.asarray([sumvals[6]-sumvals[5]])
    plt.errorbar(sumvals[5],np.max(gauss_comb_sum)*1.08,xerr=[xerrlow,xerrhigh],yerr=None,
                 fmt='o',lw=lthick, markersize=marksize,capsize=capsize,
                 markerfacecolor='green',ecolor='green',markeredgecolor='k',label='Center $\pm$ 68\%,95\% CI')

    xerrlow  = np.asarray([sumvals[5]-sumvals[3]])
    xerrhigh = np.asarray([sumvals[7]-sumvals[5]])
    plt.errorbar(sumvals[5],np.max(gauss_comb_sum)*1.13,xerr=[xerrlow,xerrhigh],yerr=None,
                 fmt='o',lw=lthick, markersize=marksize,capsize=capsize,
                 markerfacecolor='green',ecolor='green',markeredgecolor='k')

    plt.errorbar(prodvals[1],np.max(gauss_comb_prod)*prodscale*1.155,xerr=prodvals[2],yerr=None,
                 fmt='D',lw=lthick, markersize=marksize,capsize=capsize,
                 markerfacecolor='red',ecolor='red',markeredgecolor='k')
    plt.errorbar(prodvals[0],np.max(gauss_comb_prod)*prodscale*1.205,xerr=None,yerr=None,
                 fmt='*',lw=lthick, markersize=marksize,capsize=capsize,
                 markerfacecolor='red',ecolor='red',markeredgecolor='k')

    xerrlow  = np.asarray([prodvals[5]-prodvals[4]])
    xerrhigh = np.asarray([prodvals[6]-prodvals[5]])
    plt.errorbar(prodvals[5],np.max(gauss_comb_prod)*prodscale*1.055,xerr=[xerrlow,xerrhigh],yerr=None,
                 fmt='o',lw=lthick, markersize=marksize,capsize=capsize,
                 markerfacecolor='red',ecolor='red',markeredgecolor='k')

    xerrlow  = np.asarray([prodvals[5]-prodvals[3]])
    xerrhigh = np.asarray([prodvals[7]-prodvals[5]])
    plt.errorbar(prodvals[5],np.max(gauss_comb_prod)*prodscale*1.105,xerr=[xerrlow,xerrhigh],yerr=None,
                 fmt='o',lw=lthick, markersize=marksize,capsize=capsize,
                 markerfacecolor='red',ecolor='red',markeredgecolor='k')
    #-------------------------------------------------------
    if verbose: print '   plot individual Gaussian representations '
    for vv, val in enumerate(values):
        xvals, gauss = gaussreps[str(vv)]
        if vv == 0:
            plt.plot(xvals,gauss,ls='-',color='k',lw=lthick*1,label='Gauss rep.: $\mu$ $\pm$ '+str(sigrange)+'$\sigma$')
        else:
            plt.plot(xvals,gauss,ls='-',color='k',lw=lthick*1)
    plt.xlabel(axislabel, fontsize=Fsize)

    ylab = "Guassian representation of measurements "
    if skewparam != 0: ylab = 'Skewed '+ylab
    plt.ylabel(ylab, fontsize=Fsize)


    if fixxrange:
        xrange   = fixxrange

    if usexlog:
        plt.xlim([logxlow,xrange[1]])
        plt.xscale('log')
    else:
        plt.xlim(xrange)

    plt.ylim([0,np.max(gauss_comb_sum)*1.3])

    #-------------------------------------------------------
    if verbose: print '   plot legend '
    if legendlocation:
        # plt.errorbar(-20,-20,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='white', markersize=marksize,
        #              markerfacecolor='white',markeredgecolor = 'k',label='PA = '+str(PAstrings[0]))
        #
        #
        leg = plt.legend(fancybox=True, loc=legendlocation,prop={'size':Fsize-1},ncol=1,numpoints=1)#,
                         # bbox_to_anchor=(0.48, 1.))  # add the legend
        leg.get_frame().set_alpha(0.7)

    #-------------------------------------------------------
    if verbose: print '   Saving plot... '
    plt.savefig(outputname)
    plt.clf()
    plt.close('all')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if verbose: print ' - Returning the mean, median, and cetral values and uncertainties for sum and product "PDFs"'
    if verbose: print '   (  [valmed,valmean,valstd,m2sig,m1sig,center,p1sig,p2sig]  )'
    return sumvals, prodvals
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def gauss_skew(x, *p):
    """
    Return the PDF of a skewed normal distribution where alpha is the skewness paramter
    alpha = 0 represent a non-skewed regular Gaussian with mean mu and standard deviation sigma
    See https://en.wikipedia.org/wiki/Skew_normal_distribution for expression
    """
    mu, sigma, alpha = p

    t = (x-mu)/sigma
    A = 1.0 / ( sigma*np.sqrt(2.0*np.pi) )

    gauss_pdf = A * np.exp(-t**2.0 / 2.0)
    gauss_cdf = (1.0 + scipy.special.erf(alpha*t/np.sqrt(2.0))) / 2.

    gaussskew_pdf = 2.0 * gauss_pdf * gauss_cdf

    return gaussskew_pdf
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def correct4magnification_convolve(xvals,mean,sigma,alpha,mu,muerr,apply=False,verbose=True):
    """
    returning a skewed gaussian corrected for lensing magnification.
    Magnification is modeled as
    x
    x
    x
    the convolution of the measurement Gaussian with a
    non-skewed Gaussian for the magnification to account for uncetainties in magnification

    --- INPUT ---
    mean       Mean of measurement Gaussian
    sigma      Standard deviation of measurement Gaussian
    alpha      Skewness paramter of measurement Gaussian
    mu         Magnification to apply to measurement.
               Used as mean of Gaussian representation
    muerr      Approximate symmetric uncertainty on the magnification.
               Used as standard deviation of Gaussian representation
    apply      If apply is True the magnification provided will be _applied_
               to the measurement. I.e., instead of deconvolving the two signals,
               they will be convolved. In the case of Gaussians this means:
               mean_conv = mean1 + mean2; sigma_conv = sqrt(sigma1**2 + sigma2**2)
               (see http://www.tina-vision.net/docs/memos/2003-003.pdf)
    verbose

    """

    p_measurement       = mean,sigma,alpha
    gauss_measurement   = cgp.gauss_skew(xvals,*p_measurement)

    p_magnification     = mu,muerr,0.0
    xlow_mag            = mu-muerr*5.0
    xhigh_mag           = mu+muerr*5.0
    xvals_mag           = np.linspace(xlow_mag,xhigh_mag,100)
    gauss_magnification = cgp.gauss_skew(xvals_mag,*p_magnification)

    if apply:
        gauss_exvl_mag = np.convolve(gauss_measurement,gauss_magnification)
    else:
        gauss_exvl_mag = scipy.signal.deconvolve(gauss_measurement,gauss_magnification)


    pdb.set_trace()
    gauss_interp     = si.interp1d(xvals,gauss,bounds_error=False,fill_value=0.0)(xvals_comb)

    return gaussskew_pdf_exvl_mag
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def correct4magnification(xvals,Fobs,Ferrobs,mu,muerr,Nsigma=5.0,Npoints=100,verbose=True):
    """
    Returning a gaussian corrected for lensing magnification.
    Magnification is accounted for by determining the intrinsic true flux by

    Ftrue       = Fobs/mu
    Ferrtrue    = 1/mu * sqrt[ Fobs**2 * muerr**2 / mu**2 + Ferrobs**2 ]

    --- INPUT ---
    xvals      xrange for non-corrected Guassian representation of measurement
    Fobs       Mean of measurement Gaussian
    Ferrobs    Standard deviation of measurement Gaussian
    mu         Magnification to apply to measurement.
               Used to determine mean of Gaussian representation
    muerr      Symmetric uncertainty on the magnification.
               Used to determine standard deviation of Gaussian representation
    Nsigma     Number of sigmas with of gaussian to return
    verbose    Toggle verbosity

    """
    if verbose: print ' - Correcting fluxes for magnification '
    p_obs       = Fobs,Ferrobs,0.0
    gauss_obs   = cgp.gauss_skew(xvals,*p_obs)

    Ftrue       = float(Fobs)/float(mu)
    Ferrtrue    = 1/mu * np.sqrt( Fobs**2 * muerr**2 / mu**2 + Ferrobs**2 )

    p_true      = Ftrue,Ferrtrue,0.0
    xlow_true   = Ftrue-Ferrtrue*Nsigma
    xhigh_true  = Ftrue+Ferrtrue*Nsigma
    xvals_true  = np.linspace(xlow_true,xhigh_true,Npoints)
    gauss_true  = cgp.gauss_skew(xvals_true,*p_true)

    return xvals_true, gauss_true, gauss_obs
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def deconvolvetest():
    """

    """
    import scipy.stats
    xvals1 = np.linspace(0,5,1000)
    xvals2 = np.linspace(0,10,100)

    p1  =  2.5,0.2,0.0
    p2  =  4.0,1.0,0.0

    flux_init     = cgp.gauss_skew(xvals1,*p1)
    mag           = cgp.gauss_skew(xvals2,*p2)
    flux_obs      = np.convolve(flux_init, mag,'same')

    flux_init_sp  = scipy.stats.norm(2.5, 0.2).pdf(xvals1)
    mag_sp        = scipy.stats.norm(4.0, 1.0).pdf(xvals2)
    flux_obs_sp   = np.convolve(flux_init_sp, mag_sp,'same')

    flux_init_rec, residual = scipy.signal.deconvolve(flux_obs,mag)

    plt.axis([0, 10, 0, 1.25*np.max(flux_obs)])

    plt.plot(xvals2,mag,label='magnification',lw=3,c='b')
    plt.plot(xvals1,flux_init,label='intrinsic flux',lw=3,c='r')
    plt.plot(xvals1,flux_obs,label='observed flux',lw=3,c='g')

    plt.plot(flux_init_rec,label='intrinsic flux recovered',lw=3,ls='--',c='k')

    plt.plot(xvals1,flux_init_sp,ls=':',label='intrinsic flux (SciPy)',lw=2)
    plt.plot(xvals2,mag_sp,ls=':',label='magnification (SciPy)',lw=2)
    plt.plot(xvals1,flux_obs_sp,ls=':',label='observed flux (SciPy)',lw=2)

    plt.legend()
    plt.show()


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

