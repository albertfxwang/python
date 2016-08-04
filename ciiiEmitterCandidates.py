# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import commands
import sys
import glob
import MiGs
import pyfits
import numpy as np
import kbsutilities as kbs
import ciiiEmitterCandidates as cec
import macs2129_z6p8_MultiImgSource as mis
import equivalentwidth as EQW
import CDFS_MUSEvs3DHST as cm3
import matplotlib.pyplot as plt
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_candidates(verbose=True):
    """

    Get the sample(s) of candidates to investigate

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    cec.get_candidates()

    """
    catMUSE     = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1_MUSEsubcat.fits'
    cat3DHST    = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1_3DHSTinfo.fits'
    catMUSEorig = '/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits'


    if verbose: print ' ------------ OBJECT WITH CIII] FALLING IN 3D-HST G141 WAVELENGTH RANGE ------------ '
    cm3.get_candidates(zrange=[4.7,7.9],matchtol=0.5,catMUSE=catMUSE,cat3DHST=cat3DHST)

    if verbose: print ' ------------ OBJECTS WITH CIII] FALLING IN THE MUSE WAVELENGTH RANGE ------------ '
    dat = pyfits.open(catMUSEorig)[1].data
    zrange = [2.8,3.9]
    goodent = np.where((dat['REDSHIFT'] > zrange[0]) & (dat['REDSHIFT'] < zrange[1]))[0]

    if verbose: print ' - Found '+str(len(goodent))+' candidate objects in total'
    print dat['redshift'][goodent]
    for id in dat['UNIQUE_ID'][goodent]: print id

    if verbose: print ' ------------ OBJECTS WITH CIII] FALLING IN THE MUSE WAVELENGTH RANGE (3D-HST MATCH) ------------ '
    # use large matchtol to ensure all objects are listed irrespective of whether there is a close 3D-HST match
    cm3.get_candidates(zrange=zrange,matchtol=10.0,catMUSE=catMUSE,cat3DHST=cat3DHST)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def copy_spectra(outputdir='./',verbose=True,copy2Daswell=False,
                 specdir_MUSE='/Users/kschmidt/work/MUSE/emission_spectra_0.2/',
                 specdir_3DHST='/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1/'):
    """

    Get the sample(s) of candidates to investigate

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    cec.copy_spectra(outputdir='./spectra_CIIIcandidates/',copy2Daswell=True)

    """
    if os.path.isdir(outputdir):
        if verbose: print ' - Copying spectra to '+outputdir
    else:
        if verbose: print ' - The outputdir does not exists; doing nothing. \n   outputdir='+outputdir
        return

    cand3DHST           = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/candidateCIIIemitters.txt'
    candMUSE            = '/Users/kschmidt/work/MUSE/notes/160612/MUSE_CIII_potentialemitters.txt'
    candMUSE_3DHSTmatch = '/Users/kschmidt/work/MUSE/notes/160612/MUSE_CIII_potentialemitters_matchto3DHST.txt'

    dat_3DHST            = np.genfromtxt(cand3DHST          ,dtype=None,comments='#',names=True)
    dat_MUSE             = np.genfromtxt(candMUSE           ,dtype=None,comments='#',names=True)
    dat_MUSE_3DHSTmatch  = np.genfromtxt(candMUSE_3DHSTmatch,dtype=None,comments='#',names=True)


    if verbose: print ' - Copying MUSE spectra over'
    ids_MUSEall = np.append(dat_3DHST['ID_MUSE'],dat_MUSE['ID_MUSE'])
    for objid in ids_MUSEall:
        cpcmd = 'cp '+specdir_MUSE+'spectrum*'+str(objid)+'*.fits  '+outputdir
        cpout = commands.getoutput(cpcmd)
        if not cpout == '':
            print cpout

    if verbose: print ' - Reformatting and copying 3D-HST spectra over'
    ids_3DHSTall = np.append(dat_3DHST['ID_3DHST'],dat_MUSE_3DHSTmatch['ID_3DHST'])
    ids_MUSEall  = np.append(dat_3DHST['ID_MUSE'],dat_MUSE_3DHSTmatch['ID_MUSE'])
    for oo, objid_3dhst in enumerate(ids_3DHSTall):
        objid_muse = ids_MUSEall[oo]
        if verbose:
            idno    = oo+1
            infostr = '   Reformat 3D-HST 1D spectrum for '+str(objid_3dhst)+' / '+str(objid_muse)+\
                      '   ('+str(idno)+'/'+str(len(ids_3DHSTall))+')'
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        specs       = glob.glob(specdir_3DHST+'*'+str("%.5d" % objid_3dhst)+'*1D.fits')
        for spec in specs: # looping over the 1D spectra found (from multiple PAs)
            if len(spec) == 0:
                sys.exit('Something is not right; Did not find any spectrum for '+str(objid_3dhst)+
                         ' ('+str(objid_muse)+')\n'+spec)

            reformatname = 'spectrum_'+str(objid_muse)+'_'+spec.split('/')[-1].split('.')[0]+'_MiG1Dreformat.fits'

            cec.prep_3DHST_1Dspec4MiG1D(spec,outname=outputdir+reformatname,verbose=False)

            if copy2Daswell:
                cpcmd = 'cp '+spec.replace('1D.fits','2D.fits')+' '+\
                        outputdir+reformatname.replace('MiG1Dreformat.fits','2D.fits')
                cpout = commands.getoutput(cpcmd)
                if not cpout == '':
                    print cpout

    if verbose: print '\n - Done...'
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def prep_3DHST_1Dspec4MiG1D(specname,outname='DEFAULT',verbose=True):
    """

    Convert a 3DHST *1D.fits spectrum to

    --- INPUT ---
    sepcname         1D grism spectrum to reformat for MiG1D inspection.
    outname          Name (incl. paht) of file to save reformatted spectrum to.

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    newspec = cec.prep_3DHST_1Dspec4MiG1D('../candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1/goodss-31-G141_24503.1D.fits')

    """
    if outname == 'DEFAULT':
        outfile = './'+specname.split('/')[-1].split('.')[0]+'_MiG1Dreformat.fits'
    else:
        outfile = outname

    if verbose: print ' - Copying original spectrum ('+specname+')\n   to file containing output ('+outfile+')'
    cpcmd = 'cp '+specname+' '+outfile
    cpout = commands.getoutput(cpcmd)
    if not cpout == '':
        print cpout

    if verbose: print ' - Loading data and reformating:'
    hdulist = pyfits.open(outfile, mode='update')

    if verbose: print '   o Change header keywords'
    hdulist[1].header['TTYPE1']   = 'WAVE_AIR'

    hdulist[1].header['TTYPE2']   = 'FLUX'
    hdulist[1].header['TUNIT2']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE3']   = 'FLUXERR'
    hdulist[1].header['TUNIT3']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE4']   = 'CONTAM'
    hdulist[1].header['TUNIT4']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE5']   = 'TRACE'
    hdulist[1].header['TUNIT5']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE6']   = 'ETRACE'
    hdulist[1].header['TUNIT6']   = '1e-20 erg/s/cm2/A'

    hdulist[1].header['TTYPE7']   = 'SENSITIVITY'
    hdulist[1].header['TUNIT7']   = '1e-20 erg/s/cm2/A  /  E/S'

    if verbose: print '   o Adjust flux levels to be in cgs/[1e-20] units'
    hdulist[1].data['FLUX']       = hdulist[1].data['FLUX']/hdulist[1].data['SENSITIVITY']    *1.0e3
    hdulist[1].data['FLUXERR']    = hdulist[1].data['FLUXERR']/hdulist[1].data['SENSITIVITY'] *1.0e3
    hdulist[1].data['CONTAM']     = hdulist[1].data['CONTAM']/hdulist[1].data['SENSITIVITY']  *1.0e3
    hdulist[1].data['TRACE']      = hdulist[1].data['TRACE']/hdulist[1].data['SENSITIVITY']   *1.0e3
    hdulist[1].data['ETRACE']     = hdulist[1].data['ETRACE']/hdulist[1].data['SENSITIVITY']  *1.0e3

    if verbose: print '   o Correct sensititivity conversion'
    hdulist[1].data['SENSITIVITY']= 1.0e3 / hdulist[1].data['SENSITIVITY']

    if verbose: print ' - Storing correctiosn to output file'
    hdulist.flush()
    return outfile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def inspect_3DHST_2Dspec(matchcatalog='/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/candidateCIIIemitters.txt',
                         smooth=True,verbose=True):
    """

    Inspect the 2D grism spectra

    --- INPUT ---
    matchcatalog    A catalog with MUSE objects matched to 3D-HST; expects the format printed
                    to the screen by cm3.get_candidates()
    smooth          Smooth the 2D spectra in DS9 with the default Gaussian kernel?
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    cec.inspect_3DHST_2Dspec()

    """
    specdir = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1/'
    cm3.inspect_G141specs(matchcatalog,specdir=specdir,smooth=smooth,
                          ds9circlename='CIII]',verbose=verbose,oneobj=False)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def EW_from2Dspec(estimateEW=True,estimateEWlimit=True,Nsigmalimit=1,subtractcontam=True,normalizecont=True,
                           growbckhole=[1.0,1.0],ftype='sum',apertures='onebyone',mag='f125w',verbose=True,
                           filepath='/Users/kschmidt/work/MUSE/ciii_candidates/'):
    """
    Wrapper for estimating EW and EW limits of 2D grism spectra using the scripts in
    equivalenwidth.py and paper2_LBGsAbove7InGLASS.py

    Based on macs2129_z6p8_MultiImgSource.estimate_EWandEWlimits()

    --- INPUT ---
    estimateEW       Estimate EW returning the flux at the line position
    estimateEWlimit  Return a flux limit at the line position
    Nsigmalimit      If estimateEWlimit=True return the limit for this many sigma
    subtractcontam   If true the contamination extensions (CONTAM) will be subtracted before estimating
                     flux, backgrounds and EWs
    normalizecont    Offset spectrum by bck level when estimating line flux, i.e., taking line
                     flux from data array [fluxarr - np.median(flux_bckaperture)] instad of just [fluxarr]
    growbckhole      This keyword is provided to grow the whole in the background annulus aperture masks.
                     Give a fraction of size to grow the region to exclude in the background mask.
                     E.g. growbckhole=[0.05,0.45] will grow the line mask by 5% in the x direction and 45%
                     in the y-and direction and use that for the background annulus.
    ftype            The type of line flux estimate to perform.
                     'sum'      Simply summing (integrating) the pixels
                     'whtav'    Wheighted average (stack) multiplied by the aperture area (Npix) to 'integrate'
    apertures        Select what aparture file to use. Choices are: "lya", "civ", "oiii", "ciii"
    mag              Magnitude used in EW estimate
    verbose          Toggle verbosity
    filepath         Path where 'objectinfo.txt' and 'EW_lineandmaskdefinitions_*.txt' are located


    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    outputdic_f, outputdic_fl = cec.EW_from2Dspec(estimateEW=True,estimateEWlimit=False,apertures='ciii')

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fluxlimitvals = []
    outputdic_f   = {}
    outputdic_fl  = {}
    if estimateEW:
        if verbose: print ' - Will estimate EW (line fluxes)'
        fluxlimitvals.append(False)

    if estimateEWlimit:
        if verbose: print ' - Will estimate EW limits (line flux limits)'
        fluxlimitvals.append(True)

    if len(fluxlimitvals) == 0:
        if verbose: print ' - WARNING: Neither "estimateEW" nor "estimateEWlimit" set so returning None'
        return None
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading information for sources and EW masks'
    ObjInfo          = filepath+'Objectinfo.txt'
    ObjInfodata      = np.genfromtxt(ObjInfo,dtype=None, names=True, skip_header=3, comments='#')
    datadirfull_all  = ObjInfodata['datadir']
    filenames_all    = ObjInfodata['filename']
    IDs_all          = ObjInfodata['objid']

    if verbose: print '   Will use magnitude '+mag+' when estimating EWs'
    magAB_all        = ObjInfodata[mag]
    magABerr_all     = ObjInfodata[mag+'err']

    if apertures   == 'lya':
        linerestwave = 1216.0
    elif apertures == 'civ':
        linerestwave = 1549.0
    elif apertures == 'ciii':
        linerestwave = 1909.0
    else:
        sys.exit('Invalid choice of "apertures" keyword ("'+apertures+'"); ABORTING ')

    EQWsetup      = filepath+'EW_lineandmaskdefinitions_'+apertures+'.txt'
    EQWsetupdata  = np.genfromtxt(EQWsetup,dtype=None, names=True, skip_header=2, comments='#')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Estimate fluxes using estimator in equivalenwidth.EW_from2D()'
    for ff, fname in enumerate(filenames_all):
        IDs_EQW         = []
        twodfile        = fname
        objmag          = [magAB_all[ff],magABerr_all[ff]]
        objEWent        = np.where(EQWsetupdata['filename'] == twodfile)[0]

        if len(objEWent) == 0.0:
            if verbose: print ' - WARNING Did not find a match to '+fname+' in \n   '+EQWsetup
            flux_EQW, fluxerr_EQW, ewidth_EQW, ewidtherr_EQW, StoN_EQW = 99,99,99,99,0
        else:
            objlinepos      = [EQWsetupdata['linex'][objEWent][0],EQWsetupdata['liney'][objEWent][0]]
            objlinewidth    = [EQWsetupdata['linemask0'][objEWent][0],EQWsetupdata['linemask1'][objEWent][0]]
            objbckwidth     = [EQWsetupdata['bckmask0'][objEWent][0],EQWsetupdata['bckmask1'][objEWent][0]]
            objbckpos       = [EQWsetupdata['bckx'][objEWent][0],EQWsetupdata['bcky'][objEWent][0]]
            if verbose: print ' vvvvvvvvvvv estimate for line at '+str(objlinepos[0])+' A, '+\
                              str(objlinepos[1])+' arcsec vvvvvvvvvvv'

            for fluxlimitval in fluxlimitvals:

                flux_EQW, fluxerr_EQW, ewidth_EQW, ewidtherr_EQW, StoN_EQW \
                    = EQW.EW_from2D(datadirfull_all[ff]+twodfile,objlinepos,objlinewidth,objbckwidth,objmagobs=objmag,
                                    fluxlimit=fluxlimitval,verbose=verbose,weightmap='ivar',widthtype='coord',ftype=ftype,
                                    bckpos=objbckpos,subtractcontam=subtractcontam,bckcont=normalizecont,
                                    apertype_line=EQWsetupdata['linemasktype'][objEWent][0],linerestwave=linerestwave,
                                    apertype_bck=EQWsetupdata['bckmasktype'][objEWent][0],plotmasks=True,
                                    growbckhole=growbckhole)

                IDs_EQW.append(IDs_all[ff])

                if fluxlimitval:
                    outputdic_fl[fname] =    flux_EQW, fluxerr_EQW, ewidth_EQW, ewidtherr_EQW
                else:
                    outputdic_f[fname]  =    flux_EQW, fluxerr_EQW, ewidth_EQW, ewidtherr_EQW

    IDs_EQW      = np.sort(np.unique(np.asarray(IDs_EQW)))
    for fluxlimitval in fluxlimitvals:
        if fluxlimitval:
            outputdic_fl['IDs'] = IDs_EQW
        else:
            outputdic_f['IDs']  = IDs_EQW


    mis.print_EWresults(outputdic_f, outputdic_fl, splitstr='spectrum_', verbose=verbose)
    return outputdic_f, outputdic_fl
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_3DHSTmags(id3dhst,mags=['f125w','f140w','f160w'],verbose=True):
    """

    Get magnitudes of object from 3D-HST catalog

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    objmag, objmagerr, fullphot = cec.get_3DHSTmags(14831)

    """
    catalog = ' /Users/kschmidt/work/catalogs/skelton/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS'
    datPHOT = pyfits.open(catalog)[1].data
    objphot = datPHOT[np.where(datPHOT['id'] == id3dhst)[0]]
    if len(objphot) == 0:
        magresult = [-99]*len(mags)
        errresult = [-99]*len(mags)
    else:
        magresult = []
        errresult = []
        for magstr in mags:
            try:
                mag   = 25.0-2.5*np.log10(objphot['f_'+magstr])
                err   = (2.5/np.log(10)) * objphot['e_'+magstr] / objphot['f_'+magstr]
            except:
                sys.exit(' The magnitude value '+magstr+' does not appear to exist in the Skelton 3D-HST catalog. '
                                                        'Choose from:\n'+', '.join(datPHOT.columns.names))

            magresult.append(mag[0])
            errresult.append(err[0])

    if verbose:
        print '  id3dhst ',
        for mag in mags: print mag+'  '+mag+'_err',
        print '  '
        printstr = ' '+str(id3dhst)+'    '
        for mm in xrange(len(magresult)):
            printstr = printstr+str("%.2f" % magresult[mm])+' +/- '+str("%.2f" % errresult[mm])+'   '
        print printstr
    return np.asarray(magresult),np.asarray(errresult), objphot
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_MiG1Doutputsample(MiG1Doutput,CIIrange=[0,3],SiIVrange=[0,3],CIVrange=[0,3],CIIIrange=[2,3],
                          returnmag='f125w',matchtol=0.5,
                          plot2dir=False,returnall=False,verbose=True):
    """
    Get a sample of objects based on selections on a MiG1D inspection output


    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    MiG1Dfile = './MiG1D_MUSE_CIII_potentialemitters_inspection160609_final.txt'
    outputarr = cec.get_MiG1Doutputsample(MiG1Dfile,plot2dir=False,verbose=True,CIIrange=[0,3],SiIVrange=[0,3],CIVrange=[0,3],CIIIrange=[2,3])

    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading MiG1D output in '+MiG1Doutput
    MiG_dat, MiG_com, Mig_hdr = MiGs.load_MiGoutput(MiG1Doutput,verbose=False)
    IDs  = MiG_dat['ID']
    Nobj = len(IDs)
    if verbose: print '   Found '+str(Nobj)+' inspections in file '
    if len(np.unique(IDs)) != Nobj:
        if verbose: print '   WARNING: There a multiple inspections in file for the following ids: '
        print '   ',
        for objid in IDs:
            if len(np.where(IDs == objid))[0] > 1: print objid,
        print '   \n'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if returnall:
        if verbose: print ' - Returning info for all objects as returnall=True'
    else:
        if verbose: print ' - Returning info for objects based on requests (0=None; 1=tentative; 2=low-S/N; 3=high-S/N):'
        if verbose: print '     - CII : ',CIIrange
        if verbose: print '     - CIII: ',CIIIrange
        if verbose: print '     - CIV : ',CIVrange
        if verbose: print '     - SiIV: ',SiIVrange
        goodent = np.where((MiG_dat['CII']  >= CIIrange[0])  & (MiG_dat['CII']  <= CIIrange[1]) &
                           (MiG_dat['CIII'] >= CIIIrange[0]) & (MiG_dat['CIII'] <= CIIIrange[1]) &
                           (MiG_dat['CIV']  >= CIVrange[0])  & (MiG_dat['CIV']  <= CIVrange[1]) &
                           (MiG_dat['SiIV'] >= SiIVrange[0]) & (MiG_dat['SiIV']  <= SiIVrange[1]) )[0]
        Ngood = len(goodent)
        if verbose: print '   Found '+str(Ngood)+' candidates satisfying selection, ',

        if Ngood == 0:
            if verbose: print ' hence returning "None"s '
            outputarray = np.column_stack([[None], [None], [None], [None], [None], [None]])
            return outputarray
        else:
            if verbose: print ' hence continuing to collect info. '
            IDs = MiG_dat['ID'][goodent]
    IDout = IDs.tolist()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Getting info for objects '
    infofile    = '/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits'
    cmatch3dhst = '/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1_3DHSTinfo.fits'
    cmatchdat   = pyfits.open(cmatch3dhst)[1].data
    redshift    = []
    redshifterr = []
    field       = []
    mag         = []
    magerr      = []
    if verbose: print '\n#  MiG1D Object selection: CII='+str(CIIrange)+';  CIII='+str(CIIIrange)+\
                      ';  CIV='+str(CIVrange)+';  SiIV='+str(SiIVrange)
    if verbose: print '#  idMUSE   id3dhst   rmatch    '+returnmag+'  '+returnmag+'_err   field  redshift  redshift_err'
    for objid in IDout:
        objinfo = MiGs.get_objinfo(infofile,objid,'UNIQUE_ID')

        field.append(objinfo['FIELD_ID'][0])
        redshift.append(objinfo['REDSHIFT'][0])
        redshifterr.append(objinfo['REDSHIFT_ERR'][0])

        cmatch_ent = np.where(cmatchdat['id_MUSE'] == str(objid))[0]
        rmatch     = cmatchdat['r_match_arcsec'][cmatch_ent]
        if rmatch < matchtol:
            id3dhst = cmatchdat['id_3DHST'][cmatch_ent]
        else:
            id3dhst = [-999]
        objmag, objmagerr, fullphot = cec.get_3DHSTmags(id3dhst,mags=[returnmag],verbose=False)
        mag.append(objmag[0])
        magerr.append(objmagerr[0])

        if verbose:
            printstr = '   '+str(objid)+'  '+str(id3dhst[0])+'    '+str("%.5f" % rmatch[0])+'   '+\
                       str("%.2f" % objmag[0])+'  '+str("%.2f" % objmagerr[0])+'         '+\
                       str("%.2d" % objinfo['FIELD_ID'])+'     '+\
                       str("%.5f" % objinfo['REDSHIFT'])+'    '+str("%.5f" % objinfo['REDSHIFT_ERR'])
            print printstr
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot2dir:
        if verbose: print ' - Plotting final object sample. '
        cec.plot_MUSElya_forsample(IDout,redshift,voffsetlist=300.0,outputdir=plot2dir,verbose=verbose,wavetype='vac')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    outputarray = np.column_stack([IDout, field, redshift, redshifterr, mag, magerr])
    return outputarray

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_MUSElya(MUSEid,redshift,voffset=0.0,datadir='./spectra_CIIIcandidates/',outputdir='./',
                 yrangefull=[-200,300],plotSN=False,
                 showsky=False,skyspecNIR='/Users/kschmidt/work/MUSE/skytable.fits',wavetype='air',verbose=True):
    """
    Plotting the spectra (MUSE and 3D-HST) of the MUSE LAEs with zoom-ins on lines of interest.
    And 'all-in-one' plot of MiG1D produces for inspecion.

    --- INPUT ---
    MUSEid        MUSE id of object to plot
    redshift      Redshift to postion emission line markes ad
    voffset       Indicate velocity offset [km/s] wrt. to the emission line markers
    datadir       Directory containing spectra
    outputdir     Directory to save figure to
    yrangefull    Yrange of full-spectra overview
    plotSN        If true a S/N version (flux/fluxerr instead of flux) of the figure will be generated
    showsky       Show the estimated sky. For the MUSE range, the average sky estimated from the actual
                  MUSE cubes in the field the object was found in is shown. For the NIR part of the plots
                  the skyspecNIR will be used.
    skyspecNIR    Sky spectrum to plot outside the MUSE range
    wavetype      The wavelength to plot the sky in (default MUSE sky wavelength is air, but zLya is vacuum)
                  Choose between 'air' and 'vac' for air and vacuum wavelengths.
                  Will use this keyword to determine the wavelength to plot for the MUSE spectraum as well.
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    cec.plot_MUSElya(10306046,3.085,voffset=300,plotSN=False,yrangefull=[-400,1100],showsky=True,wavetype='vac')
    cec.plot_MUSElya(10306046,3.085,voffset=300,plotSN=True,yrangefull=[-3,20],showsky=True,wavetype='vac')

    cec.plot_MUSElya(11931070,4.836,voffset=300,plotSN=False,yrangefull=[-400,1100],showsky=True,wavetype='vac')
    cec.plot_MUSElya(11931070,4.836,voffset=300,plotSN=True,yrangefull=[-3,20],wavetype='vac')

    cec.plot_MUSElya(11205037,3.05787,voffset=300,plotSN=False,yrangefull=[-400,1200],showsky=True,wavetype='vac')

    cec.plot_MUSElya(11503085,3.7098,voffset=-400,plotSN=False,yrangefull=[-400,1200],showsky=True,wavetype='vac')

    cec.plot_MUSElya(10306046,3.085,voffset=300,plotSN=False,yrangefull=[-400,1200],showsky=True,wavetype='vac')
    cec.plot_MUSElya(10306046,3.085,voffset=300,plotSN=True,yrangefull=[-3,20],wavetype='vac')
    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Plotting figure for id_MUSE = '+str(MUSEid)
    redshiftplot = redshift + (voffset*(redshift+1.0) / 299792.458)
    if verbose: print ' - Will emission line markers using the redshift '+str("%.6f" % redshift)+\
                      ' (z~'+str("%.6f" % redshiftplot)+' including velocity offset of '+str(voffset)+'km/s)'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading spectra to plot from '+datadir
    file_MUSE     = glob.glob(datadir+'spectrum_'+str(MUSEid)+'.fits')
    dat_MUSE      = pyfits.open(file_MUSE[0])[1].data
    wave_MUSE     = dat_MUSE['WAVE_'+wavetype.upper()]
    flux_MUSE     = dat_MUSE['FLUX']
    ferr_MUSE     = dat_MUSE['FLUXERR']
    filllow_MUSE  = flux_MUSE-ferr_MUSE
    fillhigh_MUSE = flux_MUSE+ferr_MUSE
    wavecov_MUSE  = [np.min(wave_MUSE),np.max(wave_MUSE)]

    if plotSN:
        flux_MUSE = flux_MUSE/ferr_MUSE

    fileexist_3dhst = False
    wavecov_3dhst   = [0,0] #dummy wavelength range
    file_3dhst = glob.glob(datadir+'spectrum_'+str(MUSEid)+'_*MiG1Dreformat.fits')
    N3dhst     = len(file_3dhst)
    if N3dhst > 0:
        minwaves = []
        maxwaves = []
        for f3 in file_3dhst:
            dat_3dhst      = pyfits.open(f3)[1].data
            minwaves.append(np.min(dat_3dhst['WAVE_AIR']))
            maxwaves.append(np.max(dat_3dhst['WAVE_AIR']))

        wavecov_3dhst   = [np.min(np.asarray(minwaves)),np.max(np.asarray(maxwaves))]
        fileexist_3dhst = True
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Defining quanteties to plot; checking wavelength coverage of lines given redshift'
    llistdic      = MiGs.linelistdic(listversion='full') # loading line list for plots

    plot_CII      = False; checkwave = 1335.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_CII      = True

    plot_SiIVOIV  = False; checkwave = 1400.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_SiIVOIV  = True

    plot_CIV      = False; checkwave = 1549.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_CIV      = True

    plot_CIII     = False; checkwave = 1908.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_CIII     = True

    plot_CIIb     = False; checkwave = 2326.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_CIIb     = True

    plot_MgII     = False; checkwave = 2795.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_MgII     = True

    plot_OII      = False; checkwave = 3727.5*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_OII      = True

    plot_Hd       = False; checkwave = 4101.74*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_Hd       = True
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if showsky:
        if verbose: print ' - Getting sky spectra to plot '
        museinfo = pyfits.open('/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits')[1].data
        objent   = np.where(museinfo['UNIQUE_ID'] == str(MUSEid))[0]
        if len(objent) != 1:
            sys.exit(' Found'+str()+' matches to '+str(MUSEid)+' in MUSE info fits file... that is weird!')
        fieldno  = museinfo['FIELD_ID'][objent]

        filename      = glob.glob('/Users/kschmidt/work/MUSE/skyspectra/SKY*cdfs*-'+str("%.2d" % fieldno)+'*av.fits')
        skyMUSE       = pyfits.open(filename[0])[1].data

        if   wavetype.lower() == 'air':
            sky_wave_muse = skyMUSE['lambda']
        elif wavetype.lower() == 'vac':
            sky_wave_muse = kbs.convert_wavelength(skyMUSE['lambda'],version='air2vac')
        else:
            sys.exit('Invalid value '+wavetype+' of "wavetype"')

        sky_flux_muse = skyMUSE['data']

        skyNIR       = pyfits.open(skyspecNIR)[1].data
        sky_wave_nir = skyNIR['lam']*1e4
        sky_flux_nir = skyNIR['flux']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Plotting figure '
    figuresize_x = 13
    figuresize_y = 10

    specfigure = outputdir+'ciiiEmitterCandidatePlot_MUSEid'+str(MUSEid)+'_z'+str(redshift).replace('.','p')+\
                 '_voffset'+str(voffset).replace('.','p')+'.pdf'
    if plotSN: specfigure = specfigure.replace('.pdf','_SN.pdf')

    fig        = plt.figure(figsize=(figuresize_x,figuresize_y))
    Fsize      = 10
    LW         = 2
    plt.rc('text', usetex=True)                         # enabling LaTex rendering of text
    plt.rc('font', family='serif',size=Fsize)           # setting text font
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()

    left   = 0.06   # the left side of the subplots of the figure
    right  = 0.98   # the right side of the subplots of the figure
    bottom = 0.05   # the bottom of the subplots of the figure
    top    = 0.98   # the top of the subplots of the figure
    wspace = 0.20   # the amount of width reserved for blank space between subplots
    hspace = 0.20   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

    speccol           = 'blue'
    xlabel            = '$\lambda_\\textrm{'+wavetype.lower()+'}$ / [\AA]'
    ylabel            = '$f_\lambda / [10^{-20}$erg/s/cm$^2$/\\AA]'
    if plotSN:
        ylabel        = 'S/N'
    col_linemarker    = '#006600'
    wavescale         = 1.0
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.subplot(4, 3, 1) # Lya+NV
    windowcenter = 1229.0
    windowwidth  = 28.0
    xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

    plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
    if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

    try:
        fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
        fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
        yrange  = [fluxmin,1.1*fluxmax]
    except:
        yrange = [0,10]

    if showsky:
        if windowcenter*(redshift+1.0) < 1e4:
            sky_w   = sky_wave_muse
            sky_f   = sky_flux_muse
        else:
            sky_w   = sky_wave_nir
            sky_f   = sky_flux_nir

        skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
        skywave = sky_w[skyent]
        skylow  = np.zeros(len(skywave))
        skyflux = sky_f[skyent]
        skyhigh = skyflux
        if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

        plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    plt.xlim(xrange)
    plt.ylim(yrange)
    xrange_Lya = xrange
    yrange_Lya = yrange
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CII:
        plt.subplot(4, 3, 2) # CII
        windowcenter = 1335.
        windowwidth  = 10.
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
        if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

        try:
            fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            yrange  = [fluxmin,1.1*fluxmax]
        except:
            yrange = [0,10]

        if fileexist_3dhst:
            for f3 in file_3dhst:
                f3_dat   = pyfits.open(f3)[1].data
                wave     = f3_dat['WAVE_AIR']
                flux     = f3_dat['FLUX']-f3_dat['CONTAM']
                fluxerr  = f3_dat['FLUXERR']
                if plotSN: flux = flux/fluxerr
                contam   = f3_dat['CONTAM']
                try:
                    fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    yrange  = [fluxmin,1.1*fluxmax]
                except:
                    pass

                plt.plot(wave, flux, '-',alpha=0.8,color=speccol)

                if not plotSN:
                    plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

                    filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
                    fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
                    plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)

        if showsky:
            if windowcenter*(redshift+1.0) < 1e4:
                sky_w   = sky_wave_muse
                sky_f   = sky_flux_muse
            else:
                sky_w   = sky_wave_nir
                sky_f   = sky_flux_nir

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_CII = xrange
        yrange_CII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_SiIVOIV:
        plt.subplot(4, 3, 3) # SiIV+OIV
        windowcenter = 1397.0
        windowwidth  = 13
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
        if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

        try:
            fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            yrange  = [fluxmin,1.1*fluxmax]
        except:
            yrange = [0,10]

        if fileexist_3dhst:
            for f3 in file_3dhst:
                f3_dat   = pyfits.open(f3)[1].data
                wave     = f3_dat['WAVE_AIR']
                flux     = f3_dat['FLUX']-f3_dat['CONTAM']
                fluxerr  = f3_dat['FLUXERR']
                if plotSN: flux = flux/fluxerr
                contam   = f3_dat['CONTAM']
                try:
                    fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    yrange  = [fluxmin,1.1*fluxmax]
                except:
                    pass

                plt.plot(wave, flux, '-',alpha=0.8,color=speccol)

                if not plotSN:
                    plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

                    filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
                    fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
                    plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)

        if showsky:
            if windowcenter*(redshift+1.0) < 1e4:
                sky_w   = sky_wave_muse
                sky_f   = sky_flux_muse
            else:
                sky_w   = sky_wave_nir
                sky_f   = sky_flux_nir

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_SiIVOIV = xrange
        yrange_SiIVOIV = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CIV:
        plt.subplot(4, 3, 4) # CIV
        windowcenter = 1549.
        windowwidth  = 10.
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
        if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

        try:
            fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            yrange  = [fluxmin,1.1*fluxmax]
        except:
            yrange = [0,10]

        if fileexist_3dhst:
            for f3 in file_3dhst:
                f3_dat   = pyfits.open(f3)[1].data
                wave     = f3_dat['WAVE_AIR']
                flux     = f3_dat['FLUX']-f3_dat['CONTAM']
                fluxerr  = f3_dat['FLUXERR']
                if plotSN: flux = flux/fluxerr
                contam   = f3_dat['CONTAM']
                try:
                    fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    yrange  = [fluxmin,1.1*fluxmax]
                except:
                    pass

                plt.plot(wave, flux, '-',alpha=0.8,color=speccol)

                if not plotSN:
                    plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

                    filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
                    fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
                    plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)

        if showsky:
            if windowcenter*(redshift+1.0) < 1e4:
                sky_w   = sky_wave_muse
                sky_f   = sky_flux_muse
            else:
                sky_w   = sky_wave_nir
                sky_f   = sky_flux_nir

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_CIV = xrange
        yrange_CIV = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CIII:
        plt.subplot(4, 3, 5) # CIII]
        windowcenter = 1908.0
        windowwidth  = 10
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
        if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

        try:
            fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            yrange  = [fluxmin,1.1*fluxmax]
        except:
            yrange = [0,10]

        if fileexist_3dhst:
            for f3 in file_3dhst:
                f3_dat   = pyfits.open(f3)[1].data
                wave     = f3_dat['WAVE_AIR']
                flux     = f3_dat['FLUX']-f3_dat['CONTAM']
                fluxerr  = f3_dat['FLUXERR']
                if plotSN: flux = flux/fluxerr
                contam   = f3_dat['CONTAM']
                try:
                    fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    yrange  = [fluxmin,1.1*fluxmax]
                except:
                    pass

                plt.plot(wave, flux, '-',alpha=0.8,color=speccol)

                if not plotSN:
                    plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

                    filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
                    fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
                    plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)
        if showsky:
            if windowcenter*(redshift+1.0) < 1e4:
                sky_w   = sky_wave_muse
                sky_f   = sky_flux_muse
            else:
                sky_w   = sky_wave_nir
                sky_f   = sky_flux_nir

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_CIII = xrange
        yrange_CIII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CIIb:
        plt.subplot(4, 3, 6) # CIIb
        windowcenter = 2326.
        windowwidth  = 10
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
        if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

        try:
            fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            yrange  = [fluxmin,1.1*fluxmax]
        except:
            yrange = [0,10]

        if fileexist_3dhst:
            for f3 in file_3dhst:
                f3_dat   = pyfits.open(f3)[1].data
                wave     = f3_dat['WAVE_AIR']
                flux     = f3_dat['FLUX']-f3_dat['CONTAM']
                fluxerr  = f3_dat['FLUXERR']
                if plotSN: flux = flux/fluxerr
                contam   = f3_dat['CONTAM']

                try:
                    fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    yrange  = [fluxmin,1.1*fluxmax]
                except:
                    pass

                plt.plot(wave, flux, '-',alpha=0.8,color=speccol)

                if not plotSN:
                    plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

                    filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
                    fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
                    plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)

        if showsky:
            if windowcenter*(redshift+1.0) < 1e4:
                sky_w   = sky_wave_muse
                sky_f   = sky_flux_muse
            else:
                sky_w   = sky_wave_nir
                sky_f   = sky_flux_nir

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_CIIb = xrange
        yrange_CIIb = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_MgII:
        plt.subplot(4, 3, 7) # MgII
        windowcenter = 2795.
        windowwidth  = 10
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
        if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

        try:
            fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            yrange  = [fluxmin,1.1*fluxmax]
        except:
            yrange = [0,10]

        if fileexist_3dhst:
            for f3 in file_3dhst:
                f3_dat   = pyfits.open(f3)[1].data
                wave     = f3_dat['WAVE_AIR']
                flux     = f3_dat['FLUX']-f3_dat['CONTAM']
                fluxerr  = f3_dat['FLUXERR']
                if plotSN: flux = flux/fluxerr
                contam   = f3_dat['CONTAM']
                try:
                    fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    yrange  = [fluxmin,1.1*fluxmax]
                except:
                    pass

                plt.plot(wave, flux, '-',alpha=0.8,color=speccol)

                if not plotSN:
                    plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

                    filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
                    fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
                    plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)

        if showsky:
            if windowcenter*(redshift+1.0) < 1e4:
                sky_w   = sky_wave_muse
                sky_f   = sky_flux_muse
            else:
                sky_w   = sky_wave_nir
                sky_f   = sky_flux_nir

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_MgII = xrange
        yrange_MgII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_OII:
        plt.subplot(4, 3, 8) # [OII]
        windowcenter = 3727.5
        windowwidth  = 10
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
        if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

        try:
            fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            yrange  = [fluxmin,1.1*fluxmax]
        except:
            yrange = [0,10]

        if fileexist_3dhst:
            for f3 in file_3dhst:
                f3_dat   = pyfits.open(f3)[1].data
                wave     = f3_dat['WAVE_AIR']
                flux     = f3_dat['FLUX']-f3_dat['CONTAM']
                fluxerr  = f3_dat['FLUXERR']
                if plotSN: flux = flux/fluxerr
                contam   = f3_dat['CONTAM']
                try:
                    fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    yrange  = [fluxmin,1.1*fluxmax]
                except:
                    pass

                plt.plot(wave, flux, '-',alpha=0.8,color=speccol)

                if not plotSN:
                    plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

                    filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
                    fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
                    plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)

        if showsky:
            if windowcenter*(redshift+1.0) < 1e4:
                sky_w   = sky_wave_muse
                sky_f   = sky_flux_muse
            else:
                sky_w   = sky_wave_nir
                sky_f   = sky_flux_nir

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_OII = xrange
        yrange_OII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_Hd:
        plt.subplot(4, 3, 9) # Hd
        windowcenter = 4101.74
        windowwidth  = 10
        if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
        xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale

        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
        if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)

        try:
            fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
            yrange  = [fluxmin,1.1*fluxmax]
        except:
            yrange = [0,10]

        if fileexist_3dhst:
            for f3 in file_3dhst:
                f3_dat   = pyfits.open(f3)[1].data
                wave     = f3_dat['WAVE_AIR']
                flux     = f3_dat['FLUX']-f3_dat['CONTAM']
                fluxerr  = f3_dat['FLUXERR']
                if plotSN: flux = flux/fluxerr
                contam   = f3_dat['CONTAM']
                try:
                    fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
                    yrange  = [fluxmin,1.1*fluxmax]
                except:
                    pass

                plt.plot(wave, flux, '-',alpha=0.8,color=speccol)

                if not plotSN:
                    plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

                    filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
                    fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
                    fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
                    plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)

        if showsky:
            if windowcenter*(redshift+1.0) < 1e4:
                sky_w   = sky_wave_muse
                sky_f   = sky_flux_muse
            else:
                sky_w   = sky_wave_nir
                sky_f   = sky_flux_nir

            skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
            skywave = sky_w[skyent]
            skylow  = np.zeros(len(skywave))
            skyflux = sky_f[skyent]
            skyhigh = skyflux
            if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])

            plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
        # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        plt.xlim(xrange)
        plt.ylim(yrange)
        xrange_Hd = xrange
        yrange_Hd = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.subplot(4, 3, (10,12)) # Full specs
    plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)

    if showsky:
        skyent  = np.where((sky_wave_muse > np.min(wave_MUSE)) & (sky_wave_muse < np.max(wave_MUSE)))[0]
        skywave = sky_wave_muse[skyent]
        skylow  = np.zeros(len(skywave))
        skyhigh = sky_flux_muse[skyent]
        plt.fill_between(skywave,skylow+yrangefull[0],skyhigh+yrangefull[0],alpha=0.3,color='black')

    if fileexist_3dhst:
        for f3 in file_3dhst:
            f3_dat  = pyfits.open(f3)[1].data
            wave    = f3_dat['WAVE_AIR']
            flux    = f3_dat['FLUX']-f3_dat['CONTAM']
            fluxerr = f3_dat['FLUXERR']
            if plotSN: flux = flux/fluxerr
            contam  = f3_dat['CONTAM']
            plt.plot(wave, flux, '-',alpha=0.8,color=speccol)
            if not plotSN: plt.plot(wave, contam, '--',alpha=0.8,color=speccol)

            if showsky:
                skyent  = np.where((sky_wave_nir > np.min(wave)) & (sky_wave_nir < np.max(wave)))[0]
                skywave = sky_wave_nir[skyent]
                skylow  = np.zeros(len(skywave))
                skyflux = sky_flux_nir[skyent]
                skymax  = yrangefull[1]
                skyhigh = skyflux/np.max(skyflux)*skymax

                plt.fill_between(skywave,skylow+yrangefull[0],skyhigh+yrangefull[0],alpha=0.5,color='black')

    Dyrangefull = yrangefull[1]-yrangefull[0]
    # --- DASHED "ZOOM BOXES" ---
    plt.plot(xrange_Lya,np.zeros(2)+yrange_Lya[0],'-',color='black',lw=LW)
    plt.plot(xrange_Lya,np.zeros(2)+yrange_Lya[1],'-',color='black',lw=LW)
    plt.plot(np.zeros(2)+xrange_Lya[0],yrange_Lya,'-',color='black',lw=LW)
    plt.plot(np.zeros(2)+xrange_Lya[1],yrange_Lya,'-',color='black',lw=LW)
    plt.text(np.mean(np.asarray(xrange_Lya)),yrange_Lya[1]+0.03*Dyrangefull,'Ly$\\alpha$',
             color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_CII:
        plt.plot(xrange_CII,np.zeros(2)+yrange_CII[0],'-',color='black',lw=LW)
        plt.plot(xrange_CII,np.zeros(2)+yrange_CII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CII[0],yrange_CII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CII[1],yrange_CII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_CII)),yrange_CII[1]+0.03*Dyrangefull,'CII',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_SiIVOIV:
        plt.plot(xrange_SiIVOIV,np.zeros(2)+yrange_SiIVOIV[0],'-',color='black',lw=LW)
        plt.plot(xrange_SiIVOIV,np.zeros(2)+yrange_SiIVOIV[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_SiIVOIV[0],yrange_SiIVOIV,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_SiIVOIV[1],yrange_SiIVOIV,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_SiIVOIV)),yrange_SiIVOIV[1]+0.03*Dyrangefull,'SiIV+OIV]',rotation='vertical',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_CIV:
        plt.plot(xrange_CIV,np.zeros(2)+yrange_CIV[0],'-',color='black',lw=LW)
        plt.plot(xrange_CIV,np.zeros(2)+yrange_CIV[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIV[0],yrange_CIV,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIV[1],yrange_CIV,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_CIV)),yrange_CIV[1]+0.03*Dyrangefull,'CIV',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_CIII:
        plt.plot(xrange_CIII,np.zeros(2)+yrange_CIII[0],'-',color='black',lw=LW)
        plt.plot(xrange_CIII,np.zeros(2)+yrange_CIII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIII[0],yrange_CIII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIII[1],yrange_CIII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_CIII)),yrange_CIII[1]+0.03*Dyrangefull,'CIII]',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_CIIb:
        plt.plot(xrange_CIIb,np.zeros(2)+yrange_CIIb[0],'-',color='black',lw=LW)
        plt.plot(xrange_CIIb,np.zeros(2)+yrange_CIIb[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIIb[0],yrange_CIIb,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_CIIb[1],yrange_CIIb,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_CIIb)),yrange_CIIb[1]+0.03*Dyrangefull,'CII]',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_MgII:
        plt.plot(xrange_MgII,np.zeros(2)+yrange_MgII[0],'-',color='black',lw=LW)
        plt.plot(xrange_MgII,np.zeros(2)+yrange_MgII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_MgII[0],yrange_MgII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_MgII[1],yrange_MgII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_MgII)),yrange_MgII[1]+0.03*Dyrangefull,'MgII',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_OII:
        plt.plot(xrange_OII,np.zeros(2)+yrange_OII[0],'-',color='black',lw=LW)
        plt.plot(xrange_OII,np.zeros(2)+yrange_OII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_OII[0],yrange_OII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_OII[1],yrange_OII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_OII)),yrange_OII[1]+0.03*Dyrangefull,'[OII]',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_Hd:
        plt.plot(xrange_Hd,np.zeros(2)+yrange_Hd[0],'-',color='black',lw=LW)
        plt.plot(xrange_Hd,np.zeros(2)+yrange_Hd[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Hd[0],yrange_Hd,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Hd[1],yrange_Hd,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_Hd)),yrange_Hd[1]+0.03*Dyrangefull,'H$\delta$',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    # plt.xlim(np.min(wavecov_MUSE),np.max(wavecov_MUSE))
    plt.ylim(yrangefull)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.savefig(specfigure, dpi=300) # dpi = dot per inch for rasterized points
    plt.clf()
    plt.close('all')
    if verbose: print ' - Saved figure to ',specfigure
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW):
    """

    """
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    for ll in llistdic.keys():
        linedat      = llistdic[ll]
        linename     = linedat[0]
        if wavetype == 'vac':
            linewave = linedat[1]
        else:
            linewave = kbs.convert_wavelength(linedat[1],version='vac2air')
            pdb.set_trace()
        horalign     = linedat[2]
        lineposition = linewave*(redshift+1.0)/wavescale

        if (lineposition > xrange[0]) & (lineposition < xrange[1]):
            plt.plot(np.zeros(2)+lineposition,yrange,color=col_linemarker,alpha=0.7,linestyle='-',linewidth=LW)
            if horalign == 'right':
                xpos = lineposition-0.2*windowwidth
            elif horalign == 'left':
                xpos = lineposition+0.2*windowwidth
            else:
                xpos = lineposition

            if ll in ['oiv1','oiv2']:
                ypos = yrange[1]*0.85
            else:
                ypos = yrange[1]*0.95

            plt.text(xpos,ypos,linename,color=col_linemarker,size=Fsize,
                     rotation='horizontal',horizontalalignment=horalign,verticalalignment='top')

            if voffset != 0.0:
                zoffset  = voffset*(redshift+1.0) / 299792.458
                range    = np.sort(np.asarray([((redshift+zoffset)+1)*linewave/wavescale, (redshift +1)*linewave/wavescale]))
                lineymin = yrange[0]
                lineymax = yrange[1]
                plt.fill_between(range,np.zeros(2)+lineymin,np.zeros(2)+lineymax,alpha=0.2,color=col_linemarker)
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_MUSElya_forsample(MUSEidlist,redshiftlist,voffsetlist=0.0,outputdir='./',datadir='./spectra_CIIIcandidates/',
                           yrangefullflux=[-400,1200],yrangefullSN=[-3,30],wavetype='vac',verbose=True):
    """
    Wrapper to run cec.plot_MUSElya() for a sample of objects by providing a list of ids

    --- INPUT ---
    MUSEidlist      List of MUSE ids of objects to plot
    redshiftlist    List of Redshifts to postion emission line markes at
    voffset         Indicate velocity offset [km/s] wrt. to the emission line markers
    outputdir       Directory to save figure to
    datadir         Directory containing spectra to plot
    yrangefullflux  Yrange of full-spectra overview in flux figure
    yrangefullSN    Yrange of full-spectra overview in S/N figure
    verbose       Toggle verbosity

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    MUSEidlist = [10306046,11931070]
    zlist      = [3.085,4.836]
    cec.plot_MUSElya_forsample(MUSEidlist,zlist,voffsetlist=300.0,outputdir='./testplots160615/')

    """
    MUSEidlist   = np.asarray(MUSEidlist)
    redshiftlist = np.asarray(redshiftlist)
    Nobj         = len(MUSEidlist)
    if verbose: print ' - Plotting spectra of '+str(Nobj)+' objects in "MUSEidlist"'

    if (type(voffsetlist) == float) or (type(voffsetlist) == int):
        voffsetlist = np.zeros(Nobj) + voffsetlist
    else:
        voffsetlist = np.asarray(voffsetlist)

    for oo, objID in enumerate(MUSEidlist):
        objz = redshiftlist[oo]
        voff = voffsetlist[oo]

        if verbose:
            idno    = oo+1
            infostr = '   Plotting spectra for object '+str(objID)+' at z = '+str(objz)+' indicating voff = '+str(voff)+\
                      '   ('+str(idno)+'/'+str(Nobj)+')          '
            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        cec.plot_MUSElya(objID,objz,voffset=voff,plotSN=False,yrangefull=yrangefullflux,outputdir=outputdir,
                         showsky=True,verbose=False,wavetype=wavetype)
        cec.plot_MUSElya(objID,objz,voffset=voff,plotSN=True,yrangefull=yrangefullSN,outputdir=outputdir,
                         verbose=False,wavetype=wavetype)

    if verbose: print '\n - Done...'
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =