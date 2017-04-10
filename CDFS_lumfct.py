# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import CDFS_lumfct as cft
import collections
import kbsutilities as kbs
import numpy as np
import balff_createDataArray as bcda
import balff_utilities as butil
import pyfits
import pdb
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_fieldinfo(verbose=True):
    """
    Get Av values for MUSE CDFS fields

    --- INPUT ---


    --- EXAMPLE OF USE ---
    import CDFS_lumfct as cft
    cft.get_fieldinfo()

    """
    fieldinfo = collections.OrderedDict()

    fieldinfo['CDFS_1']  = {'name':'candels-cdfs-01', 'ra':53.062397  , 'dec':-27.80815506}
    fieldinfo['CDFS_2']  = {'name':'candels-cdfs-02', 'ra':53.06840134, 'dec':-27.82277679}
    fieldinfo['CDFS_3']  = {'name':'candels-cdfs-03', 'ra':53.07440948, 'dec':-27.83739662}
    fieldinfo['CDFS_4']  = {'name':'candels-cdfs-04', 'ra':53.08042145, 'dec':-27.85201454}
    fieldinfo['CDFS_5']  = {'name':'candels-cdfs-05', 'ra':53.08643341, 'dec':-27.86663437}
    fieldinfo['CDFS_6']  = {'name':'candels-cdfs-06', 'ra':53.07892227, 'dec':-27.80284119}
    fieldinfo['CDFS_7']  = {'name':'candels-cdfs-07', 'ra':53.08493423, 'dec':-27.81746101}
    fieldinfo['CDFS_8']  = {'name':'candels-cdfs-08', 'ra':53.09094238, 'dec':-27.83208084}
    fieldinfo['CDFS_9']  = {'name':'candels-cdfs-09', 'ra':53.09695435, 'dec':-27.84669876}
    fieldinfo['CDFS_10'] = {'name':'candels-cdfs-10', 'ra':53.10297012, 'dec':-27.86131859}
    fieldinfo['CDFS_11'] = {'name':'candels-cdfs-11', 'ra':53.09545135, 'dec':-27.79752731}
    fieldinfo['CDFS_12'] = {'name':'candels-cdfs-12', 'ra':53.1014595 , 'dec':-27.81214523}
    fieldinfo['CDFS_13'] = {'name':'candels-cdfs-13', 'ra':53.10747528, 'dec':-27.82676315}
    fieldinfo['CDFS_14'] = {'name':'candels-cdfs-14', 'ra':53.11348724, 'dec':-27.84138107}
    fieldinfo['CDFS_15'] = {'name':'candels-cdfs-15', 'ra':53.11950302, 'dec':-27.85599899}
    fieldinfo['CDFS_16'] = {'name':'candels-cdfs-16', 'ra':53.13603592, 'dec':-27.8506794 }
    fieldinfo['CDFS_17'] = {'name':'candels-cdfs-17', 'ra':53.15256882, 'dec':-27.84535599}
    fieldinfo['CDFS_18'] = {'name':'candels-cdfs-18', 'ra':53.1690979 , 'dec':-27.84003258}
    fieldinfo['CDFS_19'] = {'name':'candels-cdfs-19', 'ra':53.18562698, 'dec':-27.83470535}
    fieldinfo['CDFS_20'] = {'name':'candels-cdfs-20', 'ra':53.20215225, 'dec':-27.82937813}
    fieldinfo['CDFS_21'] = {'name':'candels-cdfs-21', 'ra':53.21867752, 'dec':-27.82404709}
    fieldinfo['CDFS_22'] = {'name':'candels-cdfs-22', 'ra':53.13002014, 'dec':-27.83606148}
    fieldinfo['CDFS_23'] = {'name':'candels-cdfs-23', 'ra':53.14654922, 'dec':-27.83073997}
    fieldinfo['CDFS_24'] = {'name':'candels-cdfs-24', 'ra':53.16307449, 'dec':-27.82541656}

    return fieldinfo
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_AVs(filter='F606W',verbose=True):
    """
    Get Av values for MUSE CDFS fields

    --- INPUT ---
    MUSEidlist

    --- EXAMPLE OF USE ---
    import CDFS_lumfct as cft
    cft.get_AVs()

    """

    fieldinfo = cft.get_fieldinfo()

    for ff in fieldinfo.keys():
        if verbose: print ' - Getting Av in filter '+filter+' for '+fieldinfo[ff]['name']
        Av,Ebv = kbs.getAv(fieldinfo[ff]['ra'],fieldinfo[ff]['dec'],filter)
        print Av, Ebv
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_objinfofits(musecat='default',photomatch='default',matchtol=0.5,photocat='default',photomag='f606w',
                       outputfitsname = './balff_data/objects_info.fits',cosmology='WMAP7BAOH0',verbose=True,
                       objecttype='Lya'):
    """
    Assemble the objects_info.fits file for BALFF

    --- INPUT ---


    --- EXAMPLE OF USE ---
    import CDFS_lumfct as cft
    cft.create_objinfofits()

    """
    if verbose: print ' - Will assemble fits table for '+objecttype+' objects'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading data'
    if musecat == 'default':
        musecat = '/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits'
    if photomatch == 'default':
        photomatch = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1_3DHSTinfo.fits'
    if photocat == 'default':
        photocat = '/Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'

    musedat       = pyfits.open(musecat)[1].data
    photomatchdat = pyfits.open(photomatch)[1].data
    photodat      = pyfits.open(photocat)[1].data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up output arrays'
    objnames        = np.array([])
    fields          = np.array([])
    Lobjs           = np.array([])
    Lobjerrs        = np.array([])
    Lfieldlims      = np.array([])
    Lfieldlimerrs   = np.array([])
    Llimsig         = 3.0 # Sigmas field limiting magnitudes correspond to
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Looping over objects to extract photometric info, calucate L and fill arrays'
    for objid in musedat['UNIQUE_ID']:
        matchent = np.where(photomatchdat['id_MUSE'] == objid)[0]
        museent  = np.where(musedat['UNIQUE_ID'] == objid)[0]
        photid   = photomatchdat['id_3DHST'][matchent]
        photent  = np.where(photodat['id'] == photid)[0]

        if (objecttype == 'all') or (musedat['IDENT_STRONGEST'][museent] == objecttype):
            Nmatch   = len(matchent)
            rmatch   = photomatchdat['r_match_arcsec'][matchent]

            if Nmatch > 1:
                matchent = matchent(np.where(rmatch == min(rmatch))[0])
                rmatch   = rmatch[matchent]
                if verbose: print '   '+objid+': More than 1 match; using match with r_match = '+str(rmatch)

            if rmatch > matchtol:
                if verbose: print '   '+objid+': Rmatch = '+str(rmatch)+' > '+str(matchtol)+' so skipping it'
                pass
            else:
                if verbose: print '   '+objid+': Extracting data...'
                # - - - MUSE CAT DATA - - -
                ra           = musedat['RA'][museent]
                dec          = musedat['DEC'][museent]
                zobj         = musedat['REDSHIFT'][museent]
                zobjerr      = musedat['REDSHIFT_ERR'][museent]
                field        = 'candels-cdfs-'+str("%.2d" % musedat['FIELD_ID'][museent])

                # - - - PHOTO CAT DATA - - -
                flux         = photodat['f_'+photomag][photent]
                fluxerr      = photodat['e_'+photomag][photent]

                Mapp         = 25.0-2.5*np.log10(flux)
                Mapperr      = (2.5/np.log(10)) * fluxerr/flux
                Mapperrfield = np.median(photodat['e_'+photomag])

                Mabs         = butil.magapp2abs(Mapp[0],zobj,ra,dec,Av=-99,band=photomag.upper(),
                                            cos=cosmology,verbose=False)
                Lobj         = butil.Mabs2L(Mabs,MUVsun=5.5)
                Lobjerr      = Lobj * np.log(10)/2.5 * Mapperr

                Lfieldlim    = butil.Mabs2L(-16.,MUVsun=5.5) # ======== HARDCORDE DUMMY VALUES ==============
                Lfieldlimerr = Lfieldlim * np.log(10)/2.5 * Mapperrfield

                # --- filling output arrays
                objnames        = np.append(objnames,objid)
                fields          = np.append(fields,field)
                Lobjs           = np.append(Lobjs,Lobj)
                Lobjerrs        = np.append(Lobjerrs,Lobjerr)
                Lfieldlims      = np.append(Lfieldlims,Lfieldlim)
                Lfieldlimerrs   = np.append(Lfieldlimerrs,Lfieldlimerr)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Generate fits file with object info for BALFF'
    output = bcda.write_fitsfile(objnames,fields,Lobjs,Lobjerrs,Lfieldlims,Lfieldlimerrs,Llimsig,
                                 outputname=outputfitsname,verbose=verbose)

    return outputfitsname
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
