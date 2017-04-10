# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# routines and wrappers to perform (force) flux measurements on the MUSE data cubes using LSDCat
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import sys
import pyfits
import commands
import numpy as np
import datetime
import collections
import astropy.wcs as wcs
import fluxmeasurementsMUSEcubes as fmm
import matplotlib.pyplot as plt
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def measure_fluxes(linecatalog, field='cdfs', field_id=15, SNthreshinit=1.0, SNthreshanal=1.0, cubeversion='_v1.0',
                   fhdu='MFS_DATA_DCBGC', ferrhdu='EFF_STAT', ffhdu='FILTERED_DATA',fferhdu='FILTERED_STAT',
                   ffsnhdu='SIGNALTONOISE',rmin=3, rmax=6, dataparentpath='/Volumes/DATABCKUP2/MUSE-Wide/data/',
                   plotfluxes=True, clobber=False,verbose=True):
    """
    Retrieve the flux measurements for an input LSDCat catalog (e.g., created with
    fmm.save_LSDCatFriendlyFitsFile) and return ...

    --- INPUT ---
    linecatalog      LSDCat catalog of (line) positions at which fluxes should be measured.
                     A mock catalog (without running LSDCat) can be generated with
                     fmm.save_LSDCatFriendlyFitsFile() to perform forced flux measurements.

    field='cdfs'
    field_id=15
    SNthreshinit=1.0
    SNthreshanal=1.0
    cubeversion      The version string used in the flux cube filenames
    fhdu             Fits extension containing flux values
    ferrhdu          Fits extension containing flux error
    ffhdu            Fits extension containing the fileterd flux
    fferhdu          Fits extension containing the fileterd flux error
    ffsnhdu          Fits extension containing the signal to noise
    rmin=3
    rmax=6
    dataparentpath
    plotfluxes       Set to True to generate flux comparison plots
    clobber          Set to True to overwrite output if it already exists
    verbose          Toggle verbosity


    --- EXAMPLE OF USE ---
    import fluxmeasurementsMUSEcubes as fmm
    linecatalog = './linecat4LSDCat_measure_161018.fits'
    fluxcatalog = fmm.measure_fluxes(linecatalog, field='cdfs', field_id=15, verbose=True)

    """
    # paths and filenames
    fieldname     = 'candels-'+field+'-'+str(field_id)
    #linecat_base=`basename ${linecatalog} .fits`

    field_path   = dataparentpath#+'/'+fieldname+'/'
    fluxcube     = field_path+'median_filtered_DATACUBE_'+fieldname+cubeversion+'.fits_effnoised_dcbgc.fits'
    SNcube       = field_path+'s2n_opt_v250_'+fieldname+'_v1.0.fits'
    filteredcube = field_path+'spec_cced_spat_cced_median_filtered_DATACUBE_'+fieldname+cubeversion+'_effnoised_32.fits'

    if verbose:
        print ' - Will measure fluxes using the following setup:'
        print """
    Field                 : %s
    Field ID              : %s
    Catalog               : %s
    fluxcube              : %s
    filtered cube         : %s
    S/N cube              : %s
    S/N thresh (init,anal): [%s,%s]
    Apertire radius range : [%s,%s]
    """ % (field,field_id,linecatalog,fluxcube,filteredcube,SNcube,SNthreshinit,SNthreshanal,rmin,rmax)

    if verbose: print ' - Putting together lsd_cat_measure.py from input '
    measure_cmd  = 'lsd_cat_measure.py ' \
                   ' --inputcat '+linecatalog+\
                   ' --thresh '+str(SNthreshinit)+\
                   ' --threshana '+str(SNthreshanal)+\
                   ' --fluxcube '+fluxcube+\
                   ' --fhdu '+fhdu+\
                   ' --ferrhdu '+ferrhdu+\
                   ' --filteredfluxcube '+filteredcube+\
                   ' --ffhdu '+ffhdu+\
                   ' --fferhdu '+fferhdu+\
                   ' --sncube '+SNcube+\
                   ' --ffsnhdu '+ffsnhdu+\
                   ' --rmin '+str(rmin)+\
                   ' --rmax '+str(rmax)

    if clobber:
        if verbose: print ' - NB: Clobber=True so adding "--clobber" to lsd_cat_measure.py command'
        measure_cmd = measure_cmd+' --clobber '

    if verbose: print ' - The command to spawn  :\n    '+measure_cmd+'\n'
    nowstr  = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    if verbose: print ' - Measuring fluxes (might take a while)'
    if verbose: print '   Started on '+nowstr
    if verbose: print '   --------------------------------- lsd_cat_measure.py output -----------------------------------'
    lsdcout = commands.getoutput(measure_cmd)
    if "DONE!!! Wrote FITS catalog" not in lsdcout:
        print ' >>>>>>>>>> WARNING: Problems with LSDCat flux measurement <<<<<<<<<<'
        print lsdcout
        print ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
        fluxcatalog = ''
    else:
        fluxcatalog = lsdcout.split('DONE!!! Wrote FITS catalog ')[-1].split(' ')[0]
        if plotfluxes:
            fmm.plot_LSDCatFluxes(fluxcatalog,verbose=verbose,zoom=False)
            fmm.plot_LSDCatFluxes(fluxcatalog,verbose=verbose,zoom=True)
    if verbose: print lsdcout
    if verbose: print '   -----------------------------------------------------------------------------------------------'
    nowstr  = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    if verbose: print '   Finished on '+nowstr

    return fluxcatalog
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_cubepixelpos(cube,ra,dec,wavelength,cubeextension='MFS_DATA_DCBGC',verbose=True):
    """
    Obtain the pixel position in a

    --- EXAMPLE OF USE ---

    import fluxmeasurementsMUSEcubes as fmm

    cube = '/Volumes/DATABCKUP3/MUSE//candels-cdfs-15/median_filtered_DATACUBE_candels-cdfs-15_v1.0.fits_effnoised.fits'
    ra_pix, dec_pix, lam_pix = fmm.get_cubepixelpos(cube,53.113493,-27.858925,6468.3,cubeextension='MEDFILTERED_DATA')  ; lam_pix=1375
    ra_pix, dec_pix, lam_pix = fmm.get_cubepixelpos(cube,53.113493,-27.858925,8506.61,cubeextension='MEDFILTERED_DATA') ; lam_pix=3006

    """
    cubehdr         = pyfits.open(cube)[cubeextension].header
    lam_pix         = fmm.get_cubepixelpos_wavelength(cubehdr,wavelength)
    ra_pix, dec_pix = fmm.get_cubepixelpos_spatial(cubehdr,ra,dec)

    if verbose:
        print '                            ra            dec           wavelength '
        print ' - Coordinates    :'+str("%15.8f" % ra)+' '+str("%15.8f" % dec)+' '+str("%15.4f" % wavelength)
        print ' - Pixel position :'+str("%15.4f" % ra_pix)+' '+str("%15.4f" % dec_pix)+' '+str("%15.4f" % lam_pix)

    return ra_pix, dec_pix, lam_pix
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_cubepixelpos_spatial(cubeheader,ra,dec,verbose=True):
    """
    Get pixel position in MUSE cube in spatial direction

    """
    del cubeheader['COMMENT'] # wcs doesn't like the comment "These data have been ZAPped!" so deleting it
    wcsobj          = wcs.WCS(cubeheader)
    lamdummy        = -99
    ra_pix, dec_pix = wcsobj.wcs_world2pix(np.array([[ra,dec,lamdummy]]),1)[0][0:2] # wcsobj.wcs_pix2world()

    return ra_pix, dec_pix
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_cubepixelpos_wavelength(cubeheader,wavelength,indexed=0,verbose=True):
    """
    Get pixel position in MUSE cube in wavelenght direction

    """
    resolution  = cubeheader['CD3_3'] # wavelength resolution in Angstrom
    waveref_val = cubeheader['CRVAL3']
    Nwaves      = cubeheader['NAXIS3']

    wavevector  = np.arange(waveref_val,waveref_val+Nwaves*resolution,resolution)

    wavediff    = np.abs(wavevector-wavelength)
    lam_pix     = np.where(wavediff == np.min(wavediff))

    if len(lam_pix) > 1:
        if verbose: print ' - WARNING: more than 1 wavelength slice matches minimum; using the lowest wavelength'
    lam_pix     = lam_pix[0]

    if indexed == 1:
        lam_pix = lam_pix + 1.0
    return lam_pix
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_LSDCatFriendlyFitsFile(outputname,lineIDs,objIDs,x_pix,y_pix,lam_pix,clobber=False,
                                radecwave=False,coordcube=None,linenames=None,verbose=True):
    """
    Save LSDCat friendly fits file which can be used as input for lsd_cat_measure.py to make
    forced flux measurements of lines or at expected line locations in the MUSE data cubes

    --- INPUT ---
    outputname  File name of output fits file to generate
    lineIDs     Unique IDs of the lines found in the objects. These IDs can potentially
                be 'appended' to the objIDs to make it easier to tie them to the parent object.
    objIDs      Unique IDs of objects to store in file. These IDs tie togehter lines in
                lineIDs list coming from the same objects.
    ra_pix      Pixel postion of line (SN peak) along the x direction in the cube
    dec_pix     Pixel postion of line (SN peak) along the y direction in the cube
    lam_pix     Pixel postion of line (SN peak) along the wavelength direction
    clobber     Set to true to overwrite existing output catalog
    radecwave   If True x_pix, y_pix and lam_pix are expected to be ra [deg], dec [deg] observed wavelength [A]
                These will then be converted to pixel positions in the cube given in coordcube
    coordcube   If x_pix, y_pix and lam_pix are given as ra, dec and wavelength (radecwave=True) provide
                the cube to be used for the coordinate (wcs) transformation from coordinates to pixel postions.
    linenames   Add list of string containing name of lines to add this to output table
    verbose     Toggle verbosity

    --- EXAMPLE OF USE ---
    import fluxmeasurementsMUSEcubes as fmm

    # >>> Example for 3 objects with 4,5 and 1 lines each <<<
    outputname = '/Users/kschmidt/work/MUSE/ciii_candidates/linecat4LSDCat_measure_161018.fits'
    lineIDs    = ['1150308501','1150308502','1150308503','1150388804']+['88801','88802','88803','88804','88805']+['99901']
    objIDs     = ['11503085']*4+[888]*5+[999]*1
    x_pix      = [120]*4+[222]*5+[333]*1
    y_pix      = [277]*4+[222]*5+[333]*1
    lam_pix    = [80,103,778,875]+[60,800,900,1200,1250]+[1111]

    fmm.save_LSDCatFriendlyFitsFile(outputname,lineIDs,objIDs,x_pix,y_pix,lam_pix,clobber=False)

    # >>> Example for multiple lines for a single cdfs 15 cube ra and dec. <<<
    linecat  = '/Users/kschmidt/work/MUSE/ciii_candidates/linecat4LSDCat_measure_cdfs15oiiemitter_161021.fits'
    obj_z    = 6474.0/3729.0-1.0
    obj_id   = '11503085'
    obj_ra   = 53.113493
    obj_dec  = -27.858925
    wave_obs = np.array([3726,3729,4340.47,9999,4861.33,4959,5007,6562.8])*(obj_z+1) # NB: Halpha outside MUSE range
    lnames   = np.asarray(['OII_1','OII_2','Hg','badline','Hb','OIII_1','OIII_2','Halpha'])
    objIDs   = np.asarray( [obj_id]*len(wave_obs) )
    lineIDs  = np.asarray( [obj_id+str("%.3d" % (nn+1)) for nn in xrange(len(wave_obs))] )
    ras      = np.asarray( [obj_ra]*len(wave_obs) )
    decs     = np.asarray( [obj_dec]*len(wave_obs) )
    cube     = '/Volumes/DATABCKUP3/MUSE//candels-cdfs-15/median_filtered_DATACUBE_candels-cdfs-15_v1.0.fits_effnoised.fits'

    fmm.save_LSDCatFriendlyFitsFile(linecat,lineIDs,objIDs,ras,decs,wave_obs,radecwave=True,coordcube=cube,linenames=lnames)

    fluxcatalog = fmm.measure_fluxes(linecat, field='cdfs', field_id=15, verbose=True, clobber=False)

    """
    datadic = {}

    if linenames == None:
        linenames = np.asarray(['None']*len(objIDs))

    if radecwave:
        if verbose: print ' - RA, Dec and wavelength provided; converting coordinates to pixel positions in cube:' \
                          '\n   '+coordcube
        # prevent changing input arrays outside this function
        ids   = objIDs.copy()
        lids  = lineIDs.copy()
        xx    = x_pix.copy()
        yy    = y_pix.copy()
        lam   = lam_pix.copy()
        lms   = linenames.copy()

        for ii in xrange(len(x_pix)):
            if (lam_pix[ii] > 4775.0) & (lam_pix[ii] < 9325.0):
                ra_pix, dec_pix, wave_pix = fmm.get_cubepixelpos(coordcube,xx[ii],yy[ii],lam[ii],verbose=False)
                xx[ii]   = ra_pix
                yy[ii]   = dec_pix
                lam[ii]  = wave_pix
            else:
                if verbose: print '   --WARNING-- Ignoring line at '+str(lam_pix[ii])+\
                                  'A as it is outside MUSE wavelength range (will not be in output table)'
                lam[ii] = -99

        datadic['I']          = np.delete(lids, np.where(lam == -99))
        datadic['ID']         = np.delete(ids,  np.where(lam == -99))
        datadic['X_PEAK_SN']  = np.delete(xx,   np.where(lam == -99))
        datadic['Y_PEAK_SN']  = np.delete(yy,   np.where(lam == -99))
        datadic['LINENAME']   = np.delete(lms,  np.where(lam == -99))
        datadic['Z_PEAK_SN']  = np.delete(lam,  np.where(lam == -99))
    else:
        datadic['I']          = np.asarray(lineIDs)
        datadic['ID']         = np.asarray(objIDs)
        datadic['X_PEAK_SN']  = np.asarray(x_pix)
        datadic['Y_PEAK_SN']  = np.asarray(y_pix)
        datadic['Z_PEAK_SN']  = np.asarray(lam_pix)
        datadic['LINENAME']   = np.asarray(linenames)

    Nlines                = len(datadic['I'])
    if verbose: print ' - Checking dimensions of input data (they should all have Nlines='+str(Nlines)+' entries)'
    lenarr = [len(datadic['ID']),len(datadic['X_PEAK_SN']),len(datadic['Y_PEAK_SN']),len(datadic['Z_PEAK_SN'])]

    if lenarr != [len(datadic['I'])]*4:
        sys.exit(' Dimensions of input data are not as expected ')
    else:
        if verbose: '   Input data looks good'
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up array for fits output \n   '+outputname
    columndefs = []
    fitsformat = ['A20','A20','D','D','D','A50']
    for kk, key in enumerate(['I','ID','X_PEAK_SN','Y_PEAK_SN','Z_PEAK_SN','LINENAME']):
        try:
            #dtype=[('I', 'S20'), ('ID', 'S20'), ('X_PEAK_SN', '>f8'), ('Y_PEAK_SN', '>f8'), ('Z_PEAK_SN', '>f8')])
            columndefs.append(pyfits.Column(name=key  , format=fitsformat[kk], array=datadic[key]))
        except:
            print ' ----ERROR---- in defining columns for fits file --> stopping with pdb.set_trace() to investigate'
            pdb.set_trace()

    cols     = pyfits.ColDefs(columndefs)
    tbhdu    = pyfits.new_table(cols)          # creating table header
    hdu      = pyfits.PrimaryHDU()             # creating primary (minimal) header
    thdulist = pyfits.HDUList([hdu, tbhdu])    # combine primary and table header to hdulist
    thdulist.writeto(outputname,clobber=clobber)  # write fits file (clobber=True overwrites excisting file)
    if verbose: print '   Wrote array to output file; done.'
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_LSDCatFluxes(fluxcatalog,margin=0.2,zoom=False,verbose=True):
    """
    Plot the content of an LSDCat flux catalog. Can for instance be generated with fmm.measure_fluxes()

    --- EXAMPLE OF USE ---
    import fluxmeasurementsMUSEcubes as fmm

    fluxcatalog = '/Users/kschmidt/work/MUSE/ciii_candidates/fluxAndEWmeasurements/11522116_linelist_fluxes.fits'
    fmm.plot_LSDCatFluxes(fluxcatalog,verbose=True)

    for fluxcat in glob.glob('/Users/kschmidt/work/MUSE/ciii_candidates/fluxAndEWmeasurements/*linelist_fluxes.fits'):
        fmm.plot_LSDCatFluxes(fluxcat,verbose=False)
        fmm.plot_LSDCatFluxes(fluxcat,verbose=True)

    """
    dat = pyfits.open(fluxcatalog)[1].data
    linecols = fmm.linecolors()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = fluxcatalog.replace('.fits','_f1kronVSf3kron.pdf')
    if zoom:
        plotname = plotname.replace('.pdf','_zoom.pdf')

    if verbose: print ' - Setting up and generating plot'
    fig = plt.figure(figsize=(6, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.98, bottom=0.10, top=0.98)
    Fsize    = 10
    lthick   = 1
    marksize = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    #plt.title(plotname.split('/')[-1].replace('_','\_'),fontsize=Fsize)

    dx      = np.abs(np.max(dat['F_KRON'])-np.min(dat['F_KRON']))
    dy      = np.abs(np.max(dat['F_3KRON'])-np.min(dat['F_3KRON']))
    zoombox = [0,2000]
    if zoom:
        xrange = zoombox
        yrange = zoombox
    else:
        xrange = [np.min(dat['F_KRON'])-dx*margin,np.max(dat['F_KRON'])+dx*margin]
        yrange = [np.min(dat['F_3KRON'])-dy*margin,np.max(dat['F_3KRON'])+dy*margin]

    for ll, line in enumerate(dat['linename']):
        for colkey in linecols.keys():
            if colkey in line.lower():
                linecol = linecols[colkey]

        plt.errorbar(dat['F_KRON'][ll],dat['F_3KRON'][ll],
                     xerr=dat['F_KRON_ERR'][ll],yerr=dat['F_3KRoN_ERR'][ll],
                     color=linecol,ls='o',lw=lthick)

        plt.text(dat['F_KRON'][ll],dat['F_3KRON'][ll],line,color=linecol,ha='left',va='bottom')

    plt.plot([np.min(xrange+yrange),np.max(xrange+yrange)],
             [np.min(xrange+yrange),np.max(xrange+yrange)],
             ls='--',color='k',lw=lthick)

    if not zoom:
        plt.plot(zoombox,np.zeros(2)+zoombox[0],ls='-',color='k',lw=lthick)
        plt.plot(zoombox,np.zeros(2)+zoombox[1],ls='-',color='k',lw=lthick)
        plt.plot(np.zeros(2)+zoombox[0],zoombox,ls='-',color='k',lw=lthick)
        plt.plot(np.zeros(2)+zoombox[1],zoombox,ls='-',color='k',lw=lthick)

    plt.xlabel(r'$F_{\rm{kron}}$ [1e-20 erg/s/cm$^2$]')
    plt.ylabel(r'$F_{3\rm{kron}}$ [1e-20 erg/s/cm$^2$]')

    plt.xlim(xrange)
    plt.ylim(yrange)

    #--------- LEGEND ---------
    # plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='white', markersize=marksize*2,
    #              markerfacecolor='white',markeredgecolor = 'k',label='Ground-based spec')
    #
    # leg = plt.legend(fancybox=True, loc='upper center',prop={'size':Fsize},ncol=1,numpoints=1)
    #                  #bbox_to_anchor=(1.25, 1.03))  # add the legend
    # leg.get_frame().set_alpha(0.7)
    #--------------------------

    if verbose: print '   Saving plot to',plotname
    plt.savefig(plotname)
    plt.clf()
    plt.close('all')
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def flux2ew(lineflux,linefluxerr,line_obswave,linerestwave,ABmag,ABmagerr,bandwavelength,
            dlambda=1.0,fluxunit=1e-20,verbose=True):
    """
    Convert flux mesurement to rest-frame EW using a broad band magnitude as continuum estimate

    --- INPUT ---
    lineflux           Line flux of line to estimate EW for
    linefluxerr        Uncertainty on line flux
    line_obswave       Observed wavelength of line in angstrom
    linerestwave       Rest-frame wavelength of line in angstrom
    ABmag              Broad band magnitude used to estimate continuum flux
    ABmagerr           Uncertainty on broad band magnitude used to estimate continuum flux
    bandwavelength     Central wavelength of broad band the AB magnitude for the continuum estimate comes from
    dlambda            Dispersion resolution to convert line fluxes from erg/s/cm2/A to erg/s/cm2 with.
                       The default is dlambda=1.0 corresponding to lineflux given in erg/s/cm2 already
    fluxunit           Units of flux; default is 1e-20
    verbose            Toggle verbosity

    --- EXAMPLE OF USE ---
    import fluxmeasurementsMUSEcubes as fmm
    lineflux        = 500.1
    linefluxerr     = 33.5
    line_obswave    = 6000.0
    linerestwave    = 1216.0
    ABmag           = 26.0
    ABmagerr        = 0.3
    bandwavelength  = 8140.0
    lineinfo = fmm.flux2ew(lineflux,linefluxerr,line_obswave,linerestwave,ABmag,ABmagerr,bandwavelength)

    """
    if verbose: print ' - Estimate continuum flux from AB magnitude [erg/s/cm2]'


    filterwave    = bandwavelength
    Fcont_mag     = 10**(-0.4*ABmag)*3*10**(-1.44) / filterwave**2 / fluxunit # continuum flux from photometry
    Fcont_mag_err = ABmagerr * Fcont_mag / (2.5/np.log(10)) # convert mag err to flux err

    if float(dlambda) != 1.0:
        if verbose: print ' - "Integrate" input line flux from [erg/s/cm2/A] to [erg/s/cm2]'
        flux         = lineflux * dlambda     # 'intergrating' units:  erg/s/cm2/A -> erg/s/cm2
        fluxerr      = linefluxerr  * dlambda
    else:
        if verbose: print ' - Input line flux in units of [erg/s/cm2] already, i.e. dlambda=1.0'
        flux         = lineflux
        fluxerr      = linefluxerr

    if verbose: print ' - Determine redshift (correction); z+1 = '
    restcorr   = line_obswave/linerestwave # = z+1   rest frame correction
    redshift   = restcorr - 1.0
    if verbose: print restcorr

    if verbose: print ' - Estimate EW and EWerr'
    ewidth     = flux      / Fcont_mag / restcorr   # [erg/s/cm2  / (erg/s/cm2/A) / (A/A)] = [A]
    ewidtherr  = np.sqrt( (fluxerr/flux)**2 + (Fcont_mag_err/Fcont_mag)**2) * np.abs(ewidth)

    if verbose: print '\n - Line at '+str(line_obswave)+'A ('+str(linerestwave)+'A restframe) at z = '+str(redshift)
    if verbose: print '   F_line  = '+str("%12.4f" % flux)     +'   +/- '+str("%12.4f" % fluxerr)+'  ['+str(fluxunit)+' erg/s/cm2]'
    if verbose: print '   ABmag   = '+str("%12.4f" % ABmag)    +'   +/- '+str("%12.4f" % ABmagerr)+'  [mag]'
    if verbose: print '   F_cont  = '+str("%12.4f" % Fcont_mag)+'   +/- '+str("%12.4f" % Fcont_mag_err)+'  ['+str(fluxunit)+' erg/s/cm2]'
    if verbose: print '   EW_rest = '+str("%12.4f" % ewidth)   +'   +/- '+str("%12.4f" % ewidtherr)+'  [A]'

    return Fcont_mag, Fcont_mag_err, ewidth, ewidtherr, redshift
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calculateEWs(fluxcatalogs,infodictionary,datadirectory='./',fluxunit=1e-20,
                 fluxcol='F_KRON',wavecol='LAMBDA_SN',clobber=False,verbose=True):
    """
    calculate the EWs for detections in LSDCat line catalog

    --- INPUT ---
    fluxcatalogs     List of fits catalogs containing line flux measrements to convert to EWs
    infodictionary   Dictionary containing inputs for fmm.flux2ew() for each of the lines in the
                     each of the flux catalogs. The format expected is:
                     infodictionary[fluxcatname] = {lineid:[linerestwave,ABmag,ABmagerr,bandwavelength]}
                     Here fluxcatname is the name of the flux catalog and lineid identifies the information
                     (not already in the fluxcatalog) needed to estimate the EW for the given line.
    datadirectory    Directory containing the flux catalogs (and the output catalogs)
    fluxunit         Unit of input fluxes
    fluxcol          Column in flux catalogs to use as estimate for line flux
    wavecol          Column in flux catalogs to use as estimate for line wavelength
    clobber          Overwrite output if it already exists
    verbose          Toggle verbosity

    --- EXAMPLE OF USE ---

    -- setup --
    import fluxmeasurementsMUSEcubes as fmm

    datadirectory   = './fluxAndEWmeasurements/'
    fluxcatalogs    = ['11503085_linelist_fluxes.fits']
    infodictionary  = {}
    infodictionary['11503085_linelist_fluxes.fits']  = {}

    #-- Lyalpha --
    f775w_skelton    = 1.705693
    f775werr_skelton = 0.033858
    linerestwave     = 1215.6737
    ABmag            = 25.0-2.5*np.log10(f775w_skelton)
    ABmagerr         = (2.5/np.log(10)) * f775werr_skelton/f775w_skelton
    bandwavelength   = 7.6485e+03

    infodictionary['11503085_linelist_fluxes.fits']['11503085003'] = [linerestwave,ABmag,ABmagerr,bandwavelength]

    #-- CIV --
    f850lp_skelton    = 1.631849
    f850lperr_skelton = 0.047236
    linerestwave      = 1548.195
    ABmag             = 25.0-2.5*np.log10(f850lp_skelton)
    ABmagerr          = (2.5/np.log(10)) * f850lperr_skelton/f850lp_skelton
    bandwavelength    = 9.1688e+03

    infodictionary['11503085_linelist_fluxes.fits']['11503085013'] = [linerestwave,ABmag,ABmagerr,bandwavelength]

    #-- Calucate EWs --
    lineinfo = fmm.calculateEWs(fluxcatalogs,infodictionary)

    """
    lineidcol = 'I' # Column of line IDs. Assumed to exist.
    Ncats     = len(fluxcatalogs)
    if verbose: print ' - Found '+str(Ncats)+' flux catalogs to estimate EWs for (in directory '+datadirectory+')'

    for cat in fluxcatalogs:
        if verbose: print ' - Loafing fluxes in '+cat+' using line ids in column '+lineidcol
        fluxdat = pyfits.open(datadirectory+cat)[1].data

        lineids_out      = []
        linerestwave_out = []
        lineflux_out     = []
        linefluxerr_out  = []
        line_obswave_out = []
        line_z_out       = []
        mag_out          = []
        magerr_out       = []
        magwave_out      = []
        fluxcont_out     = []
        fluxconterr_out  = []
        ewidth_out       = []
        ewidtherr_out    = []

        for ll, line in enumerate(fluxdat[lineidcol]):
            lineflux        = fluxdat[fluxcol][ll]
            linefluxerr     = fluxdat[fluxcol+'_ERR'][ll]
            line_obswave    = fluxdat[wavecol][ll]

            if line in infodictionary[cat].keys():
                if verbose: print ' -----------------------------------------------------------'
                if verbose: print '    - '+str(line)+':'
                lineflux        = fluxdat[fluxcol][ll]
                linefluxerr     = fluxdat[fluxcol+'_ERR'][ll]
                line_obswave    = fluxdat[wavecol][ll]
                linerestwave,ABmag,ABmagerr,bandwavelength = infodictionary[cat][line]

                lineinfo = fmm.flux2ew(lineflux,linefluxerr,line_obswave,linerestwave,ABmag,ABmagerr,bandwavelength,
                                       verbose=verbose,fluxunit=fluxunit)
                fluxcont, fluxconterr, ewidth, ewidtherr, redshift = lineinfo
                if verbose: print ' -----------------------------------------------------------'
            else:
                if verbose: print '   '+str(line)+': No setup in infodictionary; set values to -99s'
                fluxcont, fluxconterr, ewidth, ewidtherr, redshift = -99,-99,-99,-99,-99
                linerestwave,ABmag,ABmagerr,bandwavelength = -99,-99,-99,-99

            lineids_out.append(line)
            linerestwave_out.append(linerestwave)
            lineflux_out.append(lineflux)
            linefluxerr_out.append(linefluxerr)
            line_obswave_out.append(line_obswave)
            line_z_out.append(redshift)
            mag_out.append(ABmag)
            magerr_out.append(ABmagerr)
            magwave_out.append(bandwavelength)
            fluxcont_out.append(fluxcont)
            fluxconterr_out.append(fluxconterr)
            ewidth_out.append(ewidth)
            ewidtherr_out.append(ewidtherr)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print ' -----------------------------------------------------------'
        outputfile = datadirectory+'/'+cat.replace('.fits','_EWs.fits')
        if verbose: print ' - Writing results to '+outputfile
        c1  = pyfits.Column(name='line_id',       format='A20', unit='', array=np.asarray(lineids_out))
        c2  = pyfits.Column(name='line_waverest', format='D', unit='ANGSTROMS', array=np.asarray(linerestwave_out))
        c3  = pyfits.Column(name='line_flux',     format='D', unit=str(fluxunit)+' ERG/S/CM2', array=np.asarray(lineflux_out))
        c4  = pyfits.Column(name='line_flux_err', format='D', unit=str(fluxunit)+' ERG/S/CM2', array=np.asarray(linefluxerr_out))
        c5  = pyfits.Column(name='line_waveobs',  format='D', unit='ANGSTROMS', array=np.asarray(line_obswave_out))
        c6  = pyfits.Column(name='line_z',        format='D', unit='', array=np.asarray(line_z_out))
        c7  = pyfits.Column(name='mag',           format='D', unit='MAG', array=np.asarray(mag_out))
        c8  = pyfits.Column(name='mag_err',       format='D', unit='MAG', array=np.asarray(magerr_out))
        c9  = pyfits.Column(name='mag_wave',      format='D', unit='ANGSTROMS', array=np.asarray(magwave_out))
        c10 = pyfits.Column(name='cont_flux',     format='D', unit=str(fluxunit)+' ERG/S/CM2', array=np.asarray(fluxcont_out))
        c11 = pyfits.Column(name='cont_flux_err', format='D', unit=str(fluxunit)+' ERG/S/CM2', array=np.asarray(fluxconterr_out))
        c12 = pyfits.Column(name='ew',            format='D', unit='ANGSTROMS', array=np.asarray(ewidth_out))
        c13 = pyfits.Column(name='ew_err',        format='D', unit='ANGSTROMS', array=np.asarray(ewidtherr_out))

        coldefs = pyfits.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13])

        th = pyfits.new_table(coldefs) # creating default header

        # writing hdrkeys:'---KEY--',                             '----------------MAX LENGTH COMMENT-------------'
        th.header.append(('FCAT    ' , cat        ,'Flux catalog used as input for EW estimates'),end=True)
        th.header.append(('OBJID   ' , line[0:8]  ,'ID of object fluxes corrspond to'),end=True)

        head    = th.header
        tbHDU  = pyfits.new_table(coldefs, header=head)
        tbHDU.writeto(outputfile, clobber=clobber)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def linecolors():
    """
    Define line colors for plotting.

    --- EXAMPLE OF USE ---
    import fluxmeasurementsMUSEcubes as fmm
    linecolors = fmm.linecolors()

    """
    linecols = collections.OrderedDict()
    linecols['ly$\\alpha$'] = 'black'
    linecols['ly$\\beta$']  = 'gray'

    linecols['oiii']        = 'red'
    linecols['oiv']         = 'darkred'
    linecols['ovi']         = 'indianred'

    linecols['cii ']        = 'darkgreen'
    linecols['ciii]']       = 'green'
    linecols['civ']         = 'mediumspringgreen'

    linecols['nv']          = 'magenta'

    linecols['siiv']        = 'orange'

    linecols['heii']        = 'blue'

    linecols['test']        = 'yellow'
    return linecols
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =