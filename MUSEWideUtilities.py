# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                           Utilities for MUSE-Wide related stuff
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pyfits
import numpy as np
import MUSEWideUtilities as mwu
import sys
import scipy.ndimage
import pdb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def convert_wht2sigma(wht_file,wht_units='ivar',wht_ext=0,
                      smooth=False,smooth_type='gauss',smooth_param=4,
                      output_file='default',clobber=False,verbose=True):
    """
    Convert an wheight map to a sigma image.
    Add optional smoohting before saving with the smooth_* keywords.

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSEWideUtilities as mwu
    imgdir     = '/Volumes/DATABCKUP2/TDOSEextractions/tdose_cutouts/'
    wht_file   = imgdir+'acs_814w_candels-cdfs-02_wht_cut_v1.0_id9262_cutout6p0x6p0arcsec.fits'
    outputfile = mwu.convert_wht2sigma(wht_file,verbose=True,clobber=True,smooth=True,smooth_param=3)
    outputfile = mwu.convert_wht2sigma(wht_file,verbose=True,clobber=True,smooth=False)


    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading inverse variance wheight map'
    wht_data            = pyfits.open(wht_file)[wht_ext].data
    wht_data_shape      = wht_data.shape

    if wht_units == 'ivar':
        if verbose: print ' - Weight image is in units "ivar", i.e., inverser variamce. Turning map into sigmas'
        sigma_img = np.sqrt(1.0/wht_data)
    elif (wht_units == 'sigma') or (wht_units == 'stddev'):
        if verbose: print ' - Weight image already in sigma units "sigma" or "stddev". No conversion applied'
        sigma_img = wht_data
    else:
        sys.exit(' ---> Invalid choice of wht_units="'+wht_units+'" ')

    outstr = 'sigma'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if smooth:
        outstr = outstr+'_smooth'
        if smooth_type.lower() == 'gauss':
            outputimg = scipy.ndimage.filters.gaussian_filter(sigma_img, smooth_param,cval=0.0)
            smoothkey = 'gaussSig'+str(smooth_param).replace('.','p')+'pix'
            outstr    = outstr+'_'+smoothkey

        else:
            sys.exit(' ---> Invalid choice of smooth_type="'+smooth_type+'" ')
    else:
        smoothkey = 'None'
        outputimg = sigma_img

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if output_file.lower() == 'default':
        outfile = wht_file.replace('.fits','_'+outstr+'.fits')
    else:
        outfile = output_file

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Updating fits header with conversion information '
    sigma_hdr  = pyfits.open(wht_file)[wht_ext].header
    removekeys = ['FILENAME']
    for key in removekeys:
        if key in sigma_hdr.keys():
            sigma_hdr.remove(key)

    # writing hdrkeys:'---KEY--',                  '----------------MAX LENGTH COMMENT-------------'
    sigma_hdr.append(('WHTTYPE ',wht_units        ,'Assumed weight map units when converting '),end=True)
    sigma_hdr.append(('SMOOTH  ',smoothkey        ,'Information about sigma map smoothing'),end=True)
    #sigma_hdr['OBJECT'] = sigma_hdr['OBJECT'].replace('EXP','VARCUBE')

    mwu.save_img(outputimg,sigma_hdr,outfile,clobber=clobber,verbose=verbose)

    return outfile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_img(imagedata,header,filename,clobber=False,verbose=True):
    """
    Saving image to fits file

    """
    if verbose: print ' - Saving image to \n   '+filename
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if 'XTENSION' in header.keys():
        hduprim        = pyfits.PrimaryHDU()  # default HDU with default minimal header
        hducube        = pyfits.ImageHDU(imagedata,header=header)
        hdus           = [hduprim,hducube]
    else:
        hducube = pyfits.PrimaryHDU(imagedata,header=header)
        hdus           = [hducube]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hdulist = pyfits.HDUList(hdus)             # turn header into to hdulist
    hdulist.writeto(filename,clobber=clobber)  # write fits file (clobber=True overwrites excisting file)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =