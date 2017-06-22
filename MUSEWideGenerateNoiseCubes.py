# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pyfits
import numpy as np
import MUSEWideGenerateNoiseCubes as mgnc
import sys
import pdb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_variancecube(effvariancefile,exposurecubefile,ext_var=0,ext_exp='EXP',Nexpmax=4.0,
                        save_cube=True,clobber=False,verbose=True):
    """
    ---> use create_variancecube_bgrstat <---

    Generate a variance cube by scaling the effective variance by the exposure map

    NB: Generates variance cube from EFFNOISE*tables

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideGenerateNoiseCubes as mgnc
    exposurecubefile = '/Volumes/DATABCKUP2/MUSE-Wide/datacubes_dcbgc/DATACUBE_candels-cdfs-02_v1.0_dcbgc.fits'
    effvariancefile  = '/Volumes/DATABCKUP2/MUSE-Wide/effective_noise/EFFNOISE_5px_candels-cdfs-02_v1.0.fits'
    cube = mgnc.create_variancecube(effvariancefile,exposurecubefile)


    """
    if verbose: print ' - Loading effective variance and exposure map'
    var_eff             = pyfits.open(effvariancefile)[ext_var].data
    exposure_map        = pyfits.open(exposurecubefile)[ext_exp].data
    cube_shape          = exposure_map.shape

    if verbose: print ' - Checking wavelength dimensions (found '+str(len(var_eff))+' wavelength slices in effective variance)'
    if len(var_eff) == cube_shape[0]:
        if verbose: print '   The number of wavelength slices agrees between the exposure map and the effective variance'
    else:
        sys.exit(' ---> The number of wavelength slices does not agree between the exposure map and the effective variance')

    if verbose: print ' - Building cube of effective variance'
    var_eff_cube        = np.repeat( np.repeat(var_eff[:, np.newaxis],cube_shape[1],axis=1)[:,:,np.newaxis],cube_shape[2],axis=2)

    if verbose: print ' - Scaling effective variance by exposure map (masking 0s in exposure map)'
    mask_zeros          = (exposure_map == 0)
    # mask_invalid        = np.ma.masked_invalid(exposure_map).mask
    # comb_mask           = (mask_zeros | mask_invalid)

    inv_exp_map_masked  = np.ma.array(Nexpmax/exposure_map,mask=mask_zeros)
    var_eff_cube_masked = np.ma.array(var_eff_cube,mask=mask_zeros)

    var_cube_masked     = np.multiply(inv_exp_map_masked, var_eff_cube_masked)
    var_cube            = var_cube_masked.filled(fill_value=np.nan)

    if save_cube:
        varfilename     = effvariancefile.replace('.fits','_variancecube.fits')
        exposure_hdr    = pyfits.open(exposurecubefile)[ext_exp].header
        mgnc.save_cube(var_cube,exposure_hdr,varfilename,clobber=clobber,verbose=verbose)
        return varfilename
    else:
        return var_cube
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_variancecube_datapixstat(effvariancefile,exposurecubefile,ext_pixstat=1,col_var='mean',var_unitscale=1000,
                                    ext_exp='EXP',Nexpmax=4.0,save_cube=True,clobber=False,verbose=True):
    """
    ---> use create_variancecube_bgrstat <---

    Generate a variance cube by scaling the effective variance by the exposure map

    NB: Generates variance cube from datapixstat*tables

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideGenerateNoiseCubes as mgnc
    exposurecubefile = '/Volumes/DATABCKUP2/MUSE-Wide/datacubes_dcbgc/DATACUBE_candels-cdfs-02_v1.0_dcbgc.fits'
    effvariancefile  = '/Volumes/DATABCKUP2/MUSE-Wide/effective_noise/datapixstat_candels-cdfs-02.fits'
    cube = mgnc.create_variancecube_datapixstat(effvariancefile,exposurecubefile)

    """
    if verbose: print ' - Loading effective variance and exposure map'
    datapixstat_table   = pyfits.open(effvariancefile)[ext_pixstat].data
    var_eff             = datapixstat_table[col_var] * var_unitscale
    var_lambda          = datapixstat_table['wave']

    exposure_map        = pyfits.open(exposurecubefile)[ext_exp].data
    cube_shape          = exposure_map.shape

    if verbose: print ' - Checking wavelength dimensions (found '+str(len(var_lambda))+' wavelength slices in effective variance)'
    if len(var_eff) == cube_shape[0]:
        if verbose: print '   The number of wavelength slices agrees between the exposure map and the effective variance'
    else:
        sys.exit(' ---> The number of wavelength slices does not agree between the exposure map and the effective variance')

    if verbose: print ' - Building cube of effective variance'
    var_eff_cube        = np.repeat( np.repeat(var_eff[:, np.newaxis],cube_shape[1],axis=1)[:,:,np.newaxis],cube_shape[2],axis=2)

    if verbose: print ' - Scaling effective variance by exposure map (masking 0s in exposure map)'
    mask_zeros          = (exposure_map == 0)

    inv_exp_map_masked  = np.ma.array(Nexpmax/exposure_map,mask=mask_zeros)
    var_eff_cube_masked = np.ma.array(var_eff_cube,mask=mask_zeros)

    var_cube_masked     = np.multiply(inv_exp_map_masked, var_eff_cube_masked)
    var_cube            = var_cube_masked.filled(fill_value=np.nan)

    if save_cube:
        varfilename     = effvariancefile.replace('.fits','_variancecube.fits')
        exposure_hdr    = pyfits.open(exposurecubefile)[ext_exp].header
        mgnc.save_cube(var_cube,exposure_hdr,varfilename,clobber=clobber,verbose=verbose)
        return varfilename
    else:
        return var_cube
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_variancecube_bgrstat(effstatfile,exposurecubefile,ext_pixstat=1,col_sig='effsig',
                                    ext_exp='EXP',save_cube=True,clobber=False,verbose=True):
    """
    Generate a variance cube by scaling the effective standard noise by the exposure map

    NB: Generates variance cube from stat files in /store/data/musewide/bgrstat.tgz

        The STAT extension of /store/data/musewide/candels-*/DATACUBE_*_v1.0_dcbgc_effnoised.fits on arche
        contains this noise cube generated.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideGenerateNoiseCubes as mgnc
    exposurecubefile = '/Volumes/DATABCKUP2/MUSE-Wide/datacubes_dcbgc/DATACUBE_candels-cdfs-02_v1.0_dcbgc.fits'
    effstatfile      = '/Volumes/DATABCKUP2/MUSE-Wide/effective_noise/bgrstat/candels-cdfs-02/candels-cdfs-02_bgrstat.fits'
    cube             = mgnc.create_variancecube_bgrstat(effstatfile,exposurecubefile)

    """
    if verbose: print ' - Loading effective variance and exposure map'
    datapixstat_table   = pyfits.open(effstatfile)[ext_pixstat].data
    var_eff             = datapixstat_table[col_sig]**2.0
    var_lambda          = datapixstat_table['wave']

    exposure_map        = pyfits.open(exposurecubefile)[ext_exp].data
    cube_shape          = exposure_map.shape

    Nexpmax             = np.float(np.max(exposure_map[np.isfinite(exposure_map)]))

    if verbose: print ' - Checking wavelength dimensions (found '+str(len(var_lambda))+' wavelength slices in effective variance)'
    if len(var_eff) == cube_shape[0]:
        if verbose: print '   The number of wavelength slices agrees between the exposure map and the effective variance'
    else:
        sys.exit(' ---> The number of wavelength slices does not agree between the exposure map and the effective variance')

    if verbose: print ' - Building cube of effective variance'
    var_eff_cube        = np.repeat( np.repeat(var_eff[:, np.newaxis],cube_shape[1],axis=1)[:,:,np.newaxis],cube_shape[2],axis=2)

    if verbose: print ' - Scaling effective variance by exposure map (masking 0s in exposure map and using NmaxExp = '+\
                      str(Nexpmax)+')'
    mask_zeros          = (exposure_map == 0)

    inv_exp_map_masked  = np.ma.array(Nexpmax/exposure_map,mask=mask_zeros)
    var_eff_cube_masked = np.ma.array(var_eff_cube,mask=mask_zeros)

    var_cube_masked     = np.multiply(inv_exp_map_masked, var_eff_cube_masked)
    var_cube            = var_cube_masked.filled(fill_value=np.nan)

    if save_cube:
        varfilename     = effstatfile.replace('.fits','_variancecube.fits')
        exposure_hdr    = pyfits.open(exposurecubefile)[ext_exp].header
        mgnc.save_cube(var_cube,exposure_hdr,varfilename,clobber=clobber,verbose=verbose)
        return varfilename
    else:
        return var_cube


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_cube(var_cube,exposure_hdr,varfilename,clobber=False,verbose=True):
    """
    Saving variance cube to fits file

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideGenerateNoiseCubes as mgnc



    """
    if verbose: print ' - Saving variance cube to \n   '+varfilename
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Using header provided as header template'
    removekeys = ['EXTNAME','CHECKSUM','DATASUM','SCIDATA','HDUCLASS','HDUDOC','HDUVERS','HDUCLAS1']
    for key in removekeys:
        if key in exposure_hdr.keys():
            exposure_hdr.remove(key)

    # writing hdrkeys:   '---KEY--',              '----------------MAX LENGTH COMMENT-------------'
    exposure_hdr.append(('EXTNAME ','VARCUBE'    ,'Varinace cube from MUSEWideGenerateNoiseCubes()'),end=True)
    exposure_hdr['OBJECT'] = exposure_hdr['OBJECT'].replace('EXP','VARCUBE')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if 'XTENSION' in exposure_hdr.keys():
        hduprim        = pyfits.PrimaryHDU()  # default HDU with default minimal header
        hducube        = pyfits.ImageHDU(var_cube,header=exposure_hdr)
        hdus           = [hduprim,hducube]
    else:
        hducube = pyfits.PrimaryHDU(var_cube,header=exposure_hdr)
        hdus           = [hducube]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hdulist = pyfits.HDUList(hdus)       # turn header into to hdulist
    hdulist.writeto(varfilename,clobber=clobber)  # write fits file (clobber=True overwrites excisting file)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def download_effectivenoise(verbose=True):
    """
    Info...

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideGenerateNoiseCubes as mgnc
    cube = mgnc.loadcatalogs()


    """


    return None

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =