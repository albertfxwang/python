"""
#----------------------------
#   NAME
#----------------------------
# kbsutilities.py
#----------------------------
#   PURPOSE/DESCRIPTION
#----------------------------
# Collection of subroutines and scripts convenient for coding
#----------------------------
#   COMMENTS
#----------------------------
#
#----------------------------
#   EXAMPLES/USAGE
#----------------------------
# >>> import kbsutilities as kbs
#----------------------------
#   BUGS
#----------------------------
#
#----------------------------
#   REVISION HISTORY
#----------------------------
# 2013-01-04  started by K. B. Schmidt (UCSB)
#----------------------------
"""
#-------------------------------------------------------------------------------------------------------------
# importing modules
import fileinput
import sys
import math
import datetime
import numpy as np
from scipy.interpolate import interp1d
import astropysics
import astropysics.obstools
from astropysics import coords
from astropysics.constants import choose_cosmology
import types               # for testing types
import pdb
import commands
import time
import pyfits
import os
#-------------------------------------------------------------------------------------------------------------
__version__ = 1.0 
__author__ = "K. B. Schmidt"
#-------------------------------------------------------------------------------------------------------------
def test4float(str): 
    """testing whehter a string contains letters (i.e. whether it can be converted to a float)"""
    try:
        float(str)
        return True
    except ValueError:
        return False

#-------------------------------------------------------------------------------------------------------------
def replaceAll(file,searchExp,replaceExp): 
    """search and replace in file"""
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

#-------------------------------------------------------------------------------------------------------------
def arr2str(arr): 
    """Taking line of np.array or np.void and turning into string for easy passing to file"""
    strout = str(arr).replace(",",' ').replace("'",' ').replace('[',' ').replace(']',' ').replace('(',' ').replace(')',' ')+'  \n'
    return strout

#-------------------------------------------------------------------------------------------------------------
def pathAname(str): 
    """splitting string with path and name in two"""
    strsplit = str.split('/')               # splitting string
    name     = strsplit[-1]                 # saving filename (last entry of catsplit)
    slash    = '/'
    path=slash.join(strsplit[0:-1])+'/'     # putting path back together
    return [path,name]

#-------------------------------------------------------------------------------------------------------------
def DandTstr(): 
    """Creating a string with date and time on the format yymmddhhmmss"""

    now    = datetime.datetime.now()
    now    = str(now)
    now0   = now.split(' ')
    date   = now0[0].split('-')
    date1  = date[0].split('20')
    time   = now0[1].split(':')
    HHMM   = ''.join(time[0:2])
    SS     = time[2].split('.')
    HHMMSS = HHMM+str(SS[0])
    DandT = str(date1[1])+''.join(date[1:3])+HHMMSS
    return str(DandT)

#-------------------------------------------------------------------------------------------------------------
def divideWith0(numerator,denominator,badval): 
    """Replacing division with 0 in ratios with value instead of crashing"""
    return np.array([numeratorij/denominatorij if denominatorij != 0.0 else badval for numeratorij, denominatorij in zip(numerator,denominator)])

#-------------------------------------------------------------------------------------------------------------
def fluxlimit(mABobj,SNobj,SNcut=5,lam=10480,Dlam=1,Npix=1,Verbose=0): 
    """
    Calculating the limiting flux given aparent AB magnitude of an object and the obtained (average) S/N
    Default is a 5sigma limiting magnitude for the MOSFIRE Y band at 10480 Angstrom 
    (http://www2.keck.hawaii.edu/inst/mosfire/filters.html)
    To get limiting _emissionline_ flux the wavelength width and number-of-pixels values can be used. 

    ---- INPUT ----
    mABobj  The AB magnitude of the object
    SNobj   The SN of object at the selected wavelength

    ---- OPTIONAL INPUT ----
    SNcut   The SN limiting magnitude to calculate
             DEFAULT = 5 sigma
    lam     The wavelength at which the limiting magnitude flux is estimated
    Dlam    The widht of a line. Used to turn limiting flux into emission line flux
             DEFAULT = 1 (i.e. not a line)
    Npix    The number of pixels the lines (Dlam) spans over
             DEFAULT = 1 (i.e. not a line)
    Verbose Set to 1 for verbosity

    ---- EXAMPLES ----
    import kbs.utilities as kbs
    Flim     = kbs.fluxlimit(19.29,4,SNcut=3,lam=13000,Dlam=1,Npix=1,Verbose=1)
    Flimline = kbs.fluxlimit(20.28,8,SNcut=3,lam=17993,Dlam=14,Npix=1,Verbose=1)
    """

    zpAB  = 48.60
    mlim  = float(mABobj)-2.5*math.log10(float(SNcut)/float(SNobj))
    Fnu   = 10**((mlim+zpAB)/-2.5)                   # in erg/s/cm2/Hz
    Flam  =  2.99792458e+18 * Fnu / lam**2           # in erg/s/cm2/A

    Flam  = Flam * Dlam / math.sqrt(Npix)            # potentially changing Flam to limiting line flux

    if Verbose == 1: 
        print 'The ',SNcut,' sigma limiting flux was estimated to be ',Flam,' erg/s/cm^2/A'
        print 'corresponding to a limiting magnitude of ',mlim
    return Flam

#-------------------------------------------------------------------------------------------------------------
def getAv(RA,DEC,filter):
    """
    returns redening Av (and E(B-V)) for ra and dec (degrees; scalar or numpy arrays) for a given
    HST filter usning E(B-V) from the Schlegel dust maps (also returned)

       Avval, EBVval = kbs.getAv(51,219,'F125W')
    The extinction corrected apparent mag is then
       magband_corr = magband - Av_band

    Could also correct for extinction using:
        extlaw    = astropysics.obstools.CardelliExtinction(EBmV=extval, Rv=Rvval)
        magcorr   = extlaw.correctPhotometry(mag,bandwavelength)

    """
    if isinstance(RA,types.FloatType): 
        Nvals = 1
    elif isinstance(RA,types.IntType): 
        Nvals = 1
    else:
        Nvals = range(len(RA))

    if Nvals > 1:
        gall        = []
        galb        = []
        for ii in Nvals: # looping over RA and Decs and converting to galactic coordiantes
            gcoords = coords.ICRSCoordinates(RA[ii],DEC[ii]).convert(coords.GalacticCoordinates)
            gall.append(gcoords.l.degrees)
            galb.append(gcoords.b.degrees)
    else:
        gcoords = coords.ICRSCoordinates(RA,DEC).convert(coords.GalacticCoordinates)
        gall = gcoords.l.degrees
        galb = gcoords.b.degrees

    dustmaps    = '/Users/kasperborelloschmidt/work/python/progs/dust_getvalV0p1/fits/SFD_dust_4096_%s.fits'
    Ebv         = astropysics.obstools.get_SFD_dust(gall,galb,dustmaps,interpolate=True) # redening from Schlegel maps

    av_ebv = {} # ebv2Av values for HST filters from Larry
    av_ebv['F300X']  = 6.78362003559
    av_ebv['F475X']  = 3.79441819047
    av_ebv['F475W']  = 3.82839055809
    av_ebv['F606W']  = 3.01882984135
    av_ebv['F600LP'] = 2.24159324026
    av_ebv['F098M']  = 1.29502816006
    av_ebv['F105W']  = 1.18148250758
    av_ebv['F125W']  = 0.893036743585
    av_ebv['F160W']  = 0.633710427959

    try:
        av_ebv[filter] 
    except KeyError:
        sys.exit(':: kbs.getAv :: The filter '+filter+' is not accepted as input --> ABORTING')

    Av          = av_ebv[filter] * Ebv

    return Av,Ebv
#-------------------------------------------------------------------------------------------------------------
def magapp2abs(Mapp,zobj,RA,DEC,Av=-99,band='Jbradley2012',cos='WMAP7BAOH0'):
    """
    Converting apparent magnitude(s) into absolut magnitude(s)

    Av    : The extinction. If not given it's estimated from the Schlegel maps (time consuming)
            Note that RA and DEC is only used if Av is not given; otherwise they are 'dummys'
    
    band  : the band to do the calculations for. The default is to use the J band
            conversion used in Bradley et al. 2012. In this case the (extinction correted) 
            J-band magnitude is expected and MUV = MJ125 - 47.14 is returned. This 
            corresponds to 
                  Mabs = mobs - 5.0 * (np.log10(lumdist) - 1.0) + (2.5 * np.log10(1.0 + zobj))
            With k-correction (last term) assuming the source has a flat (beta = -2) SED using
            a 0.3 0.7 0.7 cosmology  
            NB! for band='Jbradley2012' zobj, RA, DEC and Av are all dummy values
    cos   : the cosmology to use, e.g.
            'WMAP7BAOH0' (Default) from http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7_bao_h0.cfm
            'WMAP7'                from http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7.cfm
    """    
    if band == 'Jbradley2012':
        Mabs          = np.array([Mapp - 47.14])
    else:
        cosmo = choose_cosmology(cos)
        Dlum          = coords.funcs.cosmo_z_to_dist(zobj, zerr=None, disttype='luminosity')*1e6 # luminosity distance in pc
        Kcorrection   = (2.5 * np.log10(1.0 + zobj)) # assumes source has flat (beta = -2) SED.  
                                                     # A bluer beta will likely give you an additional
                                                     # correction of about ~0.1 mag or so.
        if isinstance(Mapp,types.FloatType) and Av == -99: # if Av is -99, calculate it
            Av, Ebv = getAv(RA,DEC,band) 
        Mabs          = Mapp - 5*np.log10(Dlum)+5 + Kcorrection - Av # corrected absolut magnitude of objects
    return Mabs
#-------------------------------------------------------------------------------------------------------------
def magabs2app(Mabs,zobj,RA,DEC,Av=-99,band=None,cos='WMAP7BAOH0'):
    """
    Converting absolute magnitude(s) into apparent magnitude(s)

    Av    : The extinction. If not given it's estimated from the Schlegel maps (time consuming)
            Note that RA and DEC is only used if Av is not given; otherwise they are 'dummys'
    
    band  : the band to do the calculations for. The default is to use the J band
            conversion used in Bradley et al. 2012. In this case the (extinction correted) 
            J-band magnitude is expected and MJ125 = MUV + 47.14 is returned. This 
            corresponds to inverting 
                  Mabs = mobs - 5.0 * (np.log10(lumdist) - 1.0) + (2.5 * np.log10(1.0 + zobj))
            With k-correction (last term) assuming the source has a flat (beta = -2) SED using
            a 0.3 0.7 0.7 cosmology.
            NB! for band='Jbradley2012' zobj, RA, DEC and Av are all dummy values
    cos   : the cosmology to use, e.g.
            'WMAP7BAOH0' (Default) from http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7_bao_h0.cfm
            'WMAP7'                from http://lambda.gsfc.nasa.gov/product/map/dr4/params/lcdm_sz_lens_wmap7.cfm
    """
    if band == 'Jbradley2012':
        Mapp          = np.array([Mabs + 47.14])
    else:
        cosmo = choose_cosmology(cos) 
        Dlum          = coords.funcs.cosmo_z_to_dist(zobj, zerr=None, disttype='luminosity')*1e6 # luminosity distance in pc
        Kcorrection   = (2.5 * np.log10(1.0 + zobj)) # assumes source has flat (beta = -2) SED.  
                                                     # A bluer beta will likely give you an additional
                                                     # correction of about ~0.1 mag or so.
        if isinstance(Mabs,types.FloatType) and Av == -99: # if Av is -99, calculate it
            Av, Ebv = getAv(RA,DEC,band) 
        Mapp          = Mabs + 5*np.log10(Dlum) - 5 - Kcorrection + Av # corrected absolut magnitude of objects
    return Mapp
#-------------------------------------------------------------------------------------------------------------
def Mabs2L(Mabs,MUVsun=5.5):
    """
    Converting absolute magnitude(s) to luminosity in erg/s
    Using a default absolut magnitude of the sun (in UV) of 5.5 from http://www.ucolick.org/~cnaw/sun.html
    """
    Lsun        = 3.839e-11 # 1e44 erg/s
    Lobj        = 10**((MUVsun-Mabs)/2.5)*Lsun  # Luminosity in erg/s
    return Lobj
#-------------------------------------------------------------------------------------------------------------
def L2Mabs(Lobj,MUVsun=5.5):
    """
    Converting luminsoity 10^44 erg/s into absolute magnitude(s)
    Using a default absolut magnitude of the sun (in UV) of 5.5 from http://www.ucolick.org/~cnaw/sun.html
    """
    Lsun        = 3.839e-11 # 1e44 erg/s
    Mabs        = MUVsun - 2.5*np.log10(Lobj/Lsun)
    return Mabs
#-------------------------------------------------------------------------------------------------------------
def interpn(*args, **kw):
    """Interpolation on N-Dimensions 

    ai = interpn(x, y, z, ..., a, xi, yi, zi, ...)
    where the arrays x, y, z, ... define a rectangular grid
    and a.shape == (len(x), len(y), len(z), ...)

    KBS:
    Taken from http://projects.scipy.org/scipy/ticket/1727#comment:3
    An alternative is to use scipy.interpolate.LinearNDInterpolator
    but slow according to http://stackoverflow.com/questions/14119892/python-4d-linear-interpolation-on-a-rectangular-grid
    and had problems getting it to work

    -- OPTIONAL INPUT --
    method     the interpolation method to use. Options are 'linear','nearest', 'zero', 'slinear', 'quadratic', 'cubic'

    -- EAXMPLE --
    newy = kbs.interpn(oldx,oldy,newx)

    """
    method = kw.pop('method', 'linear')
    if kw:
        raise ValueError("Unknown arguments: " % kw.keys())
    nd = (len(args)-1)//2
    if len(args) != 2*nd+1:
        raise ValueError("Wrong number of arguments")
    q = args[:nd]
    qi = args[nd+1:]
    a = args[nd]
    for j in range(nd):
        a = interp1d(q[j], a, axis=j, kind=method)(qi[j])
    return a
#-------------------------------------------------------------------------------------------------------------
def simulate_schechter_distribution(alpha, L_star, L_min, N,trunmax=10):
    """ 
    Generate N samples from a Schechter distribution, which is like a gamma distribution 
    but with a negative alpha parameter and cut off on the left somewhere above zero so that
    it converges.
        
    If you pass in stupid enough parameters then it will get stuck in a loop forever, and it
    will be all your own fault.
        
    Based on algorithm in http://www.math.leidenuniv.nl/~gill/teaching/astro/stanSchechter.pdf

    KBS:-------------------------------------------------------------------------------------
          Code taken from https://gist.github.com/joezuntz/5056136 and modified.
          Schechter distribution with -1 < alpha+1 (k) < -0

          trunmax : To prevent an invinite loop trunmax gives the maximum allowed run time [s].
                    If this time is surpased any found entries are retured or an array of 0s
        -------------------------------------------------------------------------------------
    """
    output = []
    n      = 0
    Nvals  = N
    t0     = time.time()    
    while n<N:
        t1   = time.time()
        Lgam = np.random.gamma(scale=L_star, shape=alpha+2, size=N)  # drawing values from gamma dist with k+1
        Lcut = Lgam[Lgam>L_min]                                      # removing L values from gamma dist > L_min
        ucut = np.random.uniform(size=Lcut.size)                     # random values [0:1]
        Lval = Lcut[ucut<L_min/Lcut]                                 # only keeping L values where ucut < L_min/L
        output.append(Lval)                                          # append thes to output array
        n+=Lval.size                                                 # increase counter

        if (t1-t0) > trunmax:                                        # check that runtime is not too long
            Nvals = n                                                # set Nvals to found values
            if Nvals < 2.: 
                output.append(np.zeros(N))                           # if not even 2 values were found return array of 0s
                Nvals  = N                                           # updating Nvals
            n += N-n                                                 # make sure loop ends
    values = np.concatenate(output)[:Nvals]                          # generate output by reformatting
    return values
#-------------------------------------------------------------------------------------------------------------
def sex2deg(rasex,decsex):
    """ 
    converting ra and dec strings with sexagesimal values to float of degrees using skycoor

    """
    skycoorout = commands.getoutput('skycoor -d '+rasex+' '+decsex)    
    outsplit   = skycoorout.split()
    radeg      = float(outsplit[0])
    decdeg     = float(outsplit[1])
    return radeg, decdeg
#-------------------------------------------------------------------------------------------------------------
def correctmag4extinction(mag,ebv,band,Rvval=3.1,verbose=1):
    """
    Correcting a numpy array of magnitudes with given E(B-V) for galactic extionction
    Using the Cardelli et al. 1989 extintion law

    Example of usage:
       import kbsutilities as kbs
       mag  = np.array([26.75,26.25,25.88,26.36,27.31,26.97,27.32,25.77,27.07,26.91,27.64,26.67,27.04,27.16])
       ebv  = np.array([0.038,0.013,0.013,0.028,0.08,0.07,0.083,0.013,0.046,0.046,0.046,0.046,0.046,0.026])
       band = np.asarray(['F125W']*len(mag))
       magcorr = kbs.correctmag4extinction(mag,ebv,band,Rvval=3.1,verbose=1)
    """
    extval  = ebv
    extlaw  = astropysics.obstools.CardelliExtinction(EBmV=extval, Rv=Rvval)
    Nmag    = len(mag)
    wave    = np.zeros(Nmag)

    if verbose == 1: print 'Correcting the '+str(Nmag)+' magnitudes using Cardelli ext. law with Rv=',Rvval
    
    for ii in xrange(Nmag):
        bandobj = band[ii]

        if   bandobj == 'F606W': 
            wave[ii] = 5887.40
        elif bandobj == 'F098M':
            wave[ii] = 9864.10            
        elif bandobj == 'F125W':
            wave[ii] = 12486.00
        elif bandobj == 'F160W':
            wave[ii] = 15369.00
        else:
            sys.exit('The band '+str(bandobj)+' for magnitude number '+str(ii+1)+' is not a valid choice --> ABORTING')

    magcorr = extlaw.correctPhotometry(mag,wave)    
    if verbose == 1: 
        for jj in xrange(Nmag):
            print '   Corrected '+str(mag[jj])+' in '+band[jj]+' to    '+str(magcorr[jj])
    return magcorr
#-------------------------------------------------------------------------------------------------------------
def drawnbinom(n,p,size=1):
    """
    Drawing value from binomial distribution.
    Draw done using the prescription in Equation (C2) of Kelly et al. (2008)

    Note that the returned value is a float as opposed to the values returned by
        nbinomdraw  = scipy.stats.nbinom.rvs(N, p, size=1)
        nbinomdraw  = np.random.negative_binomial(N, p, size=1)
    This enables draws from distributions with very small probabilities where the
    long integers fail (returns -2147483648) in the above methods.
    
    Parameters
    ----------
        n : Number of objects in sample
        p : Probability of drawing an object from parent sample

    Returns
    -------
        A numpy array of size 'size' containing draws
    """
    t0   = time.time()    
    m    = np.zeros(size) # array to contain draws
    nlim = 1.0e6 # the limit deciding what approach to take

    if n <= nlim: # for very large n the u-vector becomes too large
        for ii in xrange(size):
            u = np.random.uniform(low=0.0, high=1.0, size=n)
            m[ii] = int(np.sum(np.log10(u)/np.log10(1-p)))    
    else:       # so in that case slit up in sub samples instead
        #print 'n > nlim = ',nlim
        nusub = np.floor(n/nlim) # sub arrays of size Nlim to draw
        nrest = n - nlim*nusub   # size of final subarray

        for ii in xrange(size):
            madd = 0.0
#            for jj in xrange(int(nusub)):
            for jj in xrange(1): # only doing it once and then multiplying by number of sub arrays
                usub  = np.random.uniform(low=0.0, high=1.0, size=nlim)
                msub  = int(np.sum(np.log10(usub)/np.log10(1-p)))    
                madd  = madd + msub * nusub
            if nrest > 0:
                urest = np.random.uniform(low=0.0, high=1.0, size=nrest)
                mrest = int(np.sum(np.log10(urest)/np.log10(1-p)))    
                madd  = madd + mrest
            
            m[ii]    = int(madd)

# looping is way too slow ...
#        sumval = 0.0
#        for ii in xrange(size):
#            for jj in xrange(int(n)):
#                u = np.random.uniform(low=0.0, high=1.0, size=1)
#                sumval = sumval + np.log10(u)/np.log10(1-p)
#            m[ii] = int(sumval)
        
    N = m + n
    #print 'Took ',time.time()-t0,' sec.'
    return N
#-------------------------------------------------------------------------------------------------------------
def appendfitstable(tab1,tab2,newtab='kbs_appendfitstable_results.fits'):
    """
    Appending 1 fits table to another.
    It is assumed that the two tables contain the same columns.
    see http://pythonhosted.org/pyfits/users_guide/users_table.html#appending-tables

    Note that columns with object IDs are also added, hence, the be aware of duplicate ids

    Parameters
    ----------
        tab1 : primariy fits table
        tab2 : fits table to append to tab1
        (should contain the same columns)

    Returns
    -------
        the name 'newtab' of the created table

    Example
    -------
    import kbsutilities as kbs
    tab1   = 'simulatedsamples/dataarraySim_pdistschechter_Ntot1000_k-0p5_Lstar0p5_LJlim0p1_Nobj17.fits'
    tab2   = 'simulatedsamples/dataarraySim_pdistschechter_Ntot2000_k-0p5_Lstar0p5_LJlim0p1_Nobj25.fits'
    newtab = 'simulatedsamples/testname.fits'
    output = kbs.appendfitstable(tab1,tab2,newtab=newtab)

    """
    t1     = pyfits.open(tab1)
    t2     = pyfits.open(tab2)

    nrows1 = t1[1].data.shape[0] # counting rows in t1
    nrows2 = t2[1].data.shape[0] # counting rows in t2

    nrows  = nrows1 + nrows2 # total number of rows in the table to be generated
    hdu    = pyfits.new_table(t1[1].columns, nrows=nrows)

    for name in t1[1].columns.names:
        hdu.data.field(name)[nrows1:]=t2[1].data.field(name)

    hdu.writeto(newtab,clobber=False)

    return newtab
#-------------------------------------------------------------------------------------------------------------
def confcontours(xpoints,ypoints,binx=200,biny=200):
    """
    Function estimating confidence contours for a given 2D distribution of points.

    @return: gridsigma, extent

    which can be plotted with for instance
    plt.contour(gridsigma.transpose(),[1,2,3],extent=extent,origin='lower',colors=['r','r','r'],label='contours',zorder=5)
    """
    from fast_kde import fast_kde # used to create confidence curves for contours
    xmin        = np.min(xpoints)
    xmax        = np.max(xpoints)
    ymin        = np.min(ypoints)
    ymax        = np.max(ypoints)
    extent      = [xmax,xmin,ymin,ymax]

    Nval        = binx*biny

    kde_grid    = fast_kde(ypoints,xpoints, gridsize=(binx,biny), weights=None,extents=[ymin,ymax,xmin,xmax])

    binarea     = (xmax-xmin)/binx * (ymax-ymin)/biny
    kde_int     = kde_grid * binarea # ~integrated value in grid
    kde_flat    = np.ravel(kde_int)
    sortindex   = np.argsort(kde_int,axis=None)[::-1]
    gridsigma   = np.zeros((binx,biny))

    sum = 0.0
    for ss in xrange(Nval):
        xx  = np.where(kde_int == kde_flat[sortindex[ss]])
        sum = sum + np.sum(kde_int[xx])
        if (sum < 0.68): gridsigma[xx] = 1.0
        if (sum > 0.68) and (sum < 0.95): gridsigma[xx] = 2.0
        if (sum > 0.95) and (sum < 0.99): gridsigma[xx] = 3.0

    return gridsigma, extent

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def createLaTeXtab(fitstable,verbose=False):
    """
    Turn binary fits table into LaTeX table

    --- INPUT ---
    fitstable   : path and name to fits catalog to turn into table

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    cat = 'dropoutselection131216_fullIR/hlsp_clash_hst_ir_macs0717_cat_reformat_flux_AB25_noF160Wcut_Bdrops_V09.fits'
    table = kbs.createLaTeXtab(cat,verbose=True)

    """
    if os.path.exists(fitstable):
        if verbose: print ' - Loading fitstable in ',fitstable
        dat   = pyfits.open(fitstable)
        datTB = dat[1].data
        cols  = datTB.columns.names
        Ncol  = len(cols)
        Nobj  = len(datTB[cols[0]])
        if verbose: print ' - Found ',Nobj,' objects to create table rows for'
        if verbose: print '   and ',Ncol,' columns of data'
    else:
        sys.exit('The provided catalog does not exist --> ABORTING')

    rows = '\n'
    for oo in xrange(Nobj):
        rowlist = [str(np.asarray(datTB[oo])[ii]) for ii in xrange(Ncol)]
        rows = rows + ' & '.join(rowlist)+'\\\\ \n'

    colheads = '\colhead{' + '} & \colhead{'.join(cols)+'}'
    colheads = colheads.replace('_','\\_')

    tablestr  = """
    \\tabletypesize{\\tiny}
    \\begin{deluxetable*}{%s}
    \\tablecolumns{%s}
    \\tablewidth{0pt}
    \\tablecaption{Title of table...}
    \\tablehead{%s}
    \startdata
    %s
    \enddata
    \\tablecomments{Comments...\\\\
    \\tnm{a}{table note...}\\\\
    }
    \label{tab:label}
    \end{deluxetable*}
    """ % (''.join(['c']*Ncol),Ncol,colheads,rows)

    return tablestr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def createLaTeXtab_CLASH(fitsCLASHcat,IDs,verbose=False,calcstuff=False):
    """
    Turn binary fits table into LaTeX table

    --- INPUT ---
    fitsCLASHcat   : path and name to fits CLASH catalog to use to turn into table
    IDs            : list of object ids to include in table

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    cat   = '/Users/kasperborelloschmidt/work/GLASS/MACS0717test/hlsp_clash_hst_ir_macs0717_cat.txt.FITS'
    table = kbs.createLaTeXtab_CLASH(cat,[859,1730],verbose=True)

    """
    if os.path.exists(fitsCLASHcat):
        if verbose: print ' - Loading CLASH catalog in ',fitsCLASHcat
        dat   = pyfits.open(fitsCLASHcat)
        datTB = dat[1].data

        Nobj    = len(IDs)
        entries = np.zeros(Nobj).astype(dtype=int)-99
        if verbose: print ' - Grabbing data for requested objects'
        for ii in xrange(Nobj):
            objent = np.where(datTB['id'] == IDs[ii])[0]
            if len(objent) == 0:
                if verbose: print '   No match to ID ',IDs[ii]
            elif len(objent) > 1:
                if verbose: print '   More than one match to ID ',IDs[ii],' --> skipping'
            else:
                entries[ii] = int(objent)

        if Nobj != len(np.where(entries != -99)[0]):
            sys.exit('Did not manage to extract data for all the '+str(Nobj)+' objects requested --> ABORTING')

        datTB = datTB[:][entries]
        cols  = datTB.columns.names
        Ncol  = len(cols)
        Nobj  = len(datTB[cols[0]])
        if verbose: print ' - Extracted ',Ncol,' colukmns for the requested ',Nobj,' objects.'
    else:
        sys.exit('The provided catalog does not exist --> ABORTING')

    if verbose: print ' - Putting data rows together.'
    headlist = ['ID','$\\alpha_\\textrm{J2000}$','$\\delta_\\textrm{J2000}$','F225W mag', 'F275W mag', 'F336W mag', 'F390W mag', 'F435W mag', 'F475W mag', 'F555W mag', 'F606W mag', 'F625W mag', 'F775W mag', 'F814W mag', 'F850lp mag', 'F105W mag', 'F110W mag', 'F125W mag', 'F140W mag', 'F160W mag', '$z_\\textrm{phot}$']

    Ncolfinal = len(headlist)
    colheads = '\colhead{' + '} & \colhead{'.join(headlist)+'}'

    fmt  = '%.2f' # format for numeric values
    fmtz = '%.1f' # format for numeric values
    rows = '\n'
    for oo in xrange(Nobj):
        rows = rows+str(int(datTB['id'][oo]))+' & '+ \
        str('%.6f' % datTB['ra'][oo])+' & '+ \
        str('%.6f' % datTB['dec'][oo])+' & ' \
        '$ '+str(fmt % datTB['f225w_mag'][oo])+' \\pm '+str(fmt % datTB['f225w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f275w_mag'][oo])+' \\pm '+str(fmt % datTB['f275w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f336w_mag'][oo])+' \\pm '+str(fmt % datTB['f336w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f390w_mag'][oo])+' \\pm '+str(fmt % datTB['f390w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f435w_mag'][oo])+' \\pm '+str(fmt % datTB['f435w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f475w_mag'][oo])+' \\pm '+str(fmt % datTB['f475w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f555w_mag'][oo])+' \\pm '+str(fmt % datTB['f555w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f606w_mag'][oo])+' \\pm '+str(fmt % datTB['f606w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f625w_mag'][oo])+' \\pm '+str(fmt % datTB['f625w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f775w_mag'][oo])+' \\pm '+str(fmt % datTB['f775w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f814w_mag'][oo])+' \\pm '+str(fmt % datTB['f814w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f850lp_mag'][oo])+' \\pm '+str(fmt % datTB['f850lp_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f105w_mag'][oo])+' \\pm '+str(fmt % datTB['f105w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f110w_mag'][oo])+' \\pm '+str(fmt % datTB['f110w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f125w_mag'][oo])+' \\pm '+str(fmt % datTB['f125w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f140w_mag'][oo])+' \\pm '+str(fmt % datTB['f140w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmt % datTB['f160w_mag'][oo])+' \\pm '+str(fmt % datTB['f160w_magerr'][oo])+' $ & ' + \
        '$ '+str(fmtz % datTB['zb'][oo]) + \
        ' ^{+'+str(fmtz % (datTB['zbmax'][oo]-datTB['zb'][oo])) + \
        '}_{-'+str(fmtz % (datTB['zb'][oo]-datTB['zbmin'][oo]))+' }$ \\\\ \n'

    tablestr  = """
    \\tabletypesize{\\tiny}
    \\begin{deluxetable*}{%s}
    \\tablecolumns{%s}
    \\tablewidth{0pt}
    \\tablecaption{Title of table...}
    \\tablehead{%s}
    \startdata
    %s
    \enddata
    \\tablecomments{Comments...\\\\
    \\tnm{a}{table note...}\\\\
    }
    \label{tab:label}
    \end{deluxetable*}
    """ % (''.join(['c']*Ncolfinal),Ncolfinal,colheads,rows)

    # ===== Calculate stuff for objects =====
    if calcstuff==True:
        if verbose: print ' - Calculating stuff for objects:'
        for oo in xrange(Nobj):
            mags  = np.array([datTB['f125w_mag'][oo],datTB['f140w_mag'][oo],datTB['f160w_mag'][oo]])
            errs  = np.array([datTB['f125w_magerr'][oo],datTB['f140w_magerr'][oo],datTB['f160w_magerr'][oo]])

            magav = np.average(mags[mags != -99],weights=1/errs[mags != -99]**2) # wighted average of magnitudes
            #if verbose: print '   obj,magav,zphot,mags[mags!=-99] = ',IDs[oo],magav,datTB['zb'][oo],mags[mags != -99]

            #if verbose: print IDs[oo],datTB['ra'][oo],datTB['dec'][oo]
            #if verbose: print IDs[oo],'   ',str(fmt % datTB['f105w_mag'][oo]),',',str(fmt % datTB['f105w_magerr'][oo])
            if verbose: print IDs[oo],'   ',str(fmt % datTB['f140w_mag'][oo]),',',str(fmt % datTB['f140w_magerr'][oo])


    return tablestr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def createSCImCONTAM(spec2Dfits,verbose=False,clobber=False):
    """
    Read a 2D fits spectrum from the 3D-HST/GLASS reduction and produce a
    new fits image containing the contamination subtracted 2D spectrum.

    --- INPUT ---
    spec2Dfits   : path and anme to 2D spectrum (ssumes it contains the extensions SCI and CONTAM)

    --- EXAMPLE OF USE ---
    import kbsutilities as kbs
    outputimg = kbs.createSCImCONTAM(spec2Dfits,verbose=True)

    """
    if os.path.isfile(spec2Dfits):
        if verbose: print ' - Loading fits file'
        twod       = pyfits.open(spec2Dfits)
        SCI        = twod['SCI'].data
        CONTAM     = twod['CONTAM'].data
        SCImCONTAM = SCI - CONTAM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def void2file(npvoid,outputfile,verbose=False,clobber=False):
    """
    Wirting a numpy void (from np.genfromtxt) to a text file
    """

    if os.path.isfile(outputfile) and (clobber == False):
        sys.exit(' clobber == false so not creating '+outputfile+' as it already exists')
    else:
        f = open(outputfile,'w')
        colnames = npvoid.dtype.names
        hdr = '# '+str(colnames).replace(',',' ').replace(')',' ').replace('(',' ').replace("'",' ').replace("'",' ').replace(']',' ').replace('[',' ')+'\n'
        f.write(hdr)
        for ii in xrange(len(npvoid[colnames[0]])):
            outstr = str(npvoid[ii].tolist()).replace(',',' ').replace(')',' ').replace('(',' ').replace("'",' ').replace("'",' ').replace(']',' ').replace('[',' ')
            f.write("%s\n" % outstr)

        f.close()
#-------------------------------------------------------------------------------------------------------------
#                                                  END
#-------------------------------------------------------------------------------------------------------------
