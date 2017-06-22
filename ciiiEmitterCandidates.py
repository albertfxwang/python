# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import commands
import sys
import glob
import MiGs
import pyfits
import datetime
import numpy as np
import kbsutilities as kbs
import ciiiEmitterCandidates as cec
import macs2129_z6p8_MultiImgSource as mis
import fluxmeasurementsMUSEcubes as fmm
import rxj2248_BooneBalestraSource as bbs
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

    if verbose: print ' ------------ OBJECTS WITH CIV FALLING IN THE MUSE WAVELENGTH RANGE BUT WITH CIII] OUTSIDE MUSE and 3D-HST (3D-HST MATCH) ------------ '
    zrange = [3.9,4.7]
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

    cand3DHST           = '/Users/kschmidt/work/MUSE/ciii_candidates/CIIIemitter_candidates_within3DHST.txt'
    candMUSE            = '/Users/kschmidt/work/MUSE/ciii_candidates/CIIIemitter_candidates_withinMUSE.txt'
    candMUSE_CIV        = '/Users/kschmidt/work/MUSE/ciii_candidates/CIVemitter_candidates_withinMUSE.txt'

    dat_3DHST            = np.genfromtxt(cand3DHST   ,skip_header=4,names=True,dtype=None,comments='#')
    dat_MUSE             = np.genfromtxt(candMUSE    ,skip_header=4,names=True,dtype=None,comments='#')
    dat_MUSE_CIV         = np.genfromtxt(candMUSE_CIV,skip_header=4,names=True,dtype=None,comments='#')

    ids_MUSE_all = np.append(dat_3DHST['ID_MUSE'],dat_MUSE['ID_MUSE'])
    ids_MUSE_all = np.append(ids_MUSE_all,dat_MUSE_CIV['ID_MUSE'])

    ids_3DHST_3dhstmatch = np.append(dat_3DHST['ID_3DHST'],dat_MUSE['ID_3DHST'][dat_MUSE['R_MATCH'] < 0.5])
    ids_3DHST_3dhstmatch = np.append(ids_3DHST_3dhstmatch,dat_MUSE_CIV['ID_3DHST'][dat_MUSE_CIV['R_MATCH'] < 0.5])

    ids_MUSE_3dhstmatch  = np.append(dat_3DHST['ID_MUSE'],dat_MUSE['ID_MUSE'][dat_MUSE['R_MATCH'] < 0.5])
    ids_MUSE_3dhstmatch  = np.append(ids_MUSE_3dhstmatch,dat_MUSE_CIV['ID_MUSE'][dat_MUSE_CIV['R_MATCH'] < 0.5])

    if verbose: print ' - Copying MUSE spectra over'
    for objid in ids_MUSE_all:
        cpcmd = 'cp '+specdir_MUSE+'spectrum*'+str(objid)+'*.fits  '+outputdir
        cpout = commands.getoutput(cpcmd)
        if not cpout == '':
            print cpout

    if verbose: print ' - Reformatting and copying 3D-HST spectra over'
    for oo, objid_3dhst in enumerate(ids_3DHST_3dhstmatch):
        objid_muse = ids_MUSE_3dhstmatch[oo]
        if verbose:
            idno    = oo+1
            infostr = '   Reformat 3D-HST 1D spectrum for '+str(objid_3dhst)+' / '+str(objid_muse)+\
                      '   ('+str(idno)+'/'+str(len(ids_MUSE_3dhstmatch))+')'
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
def get_LSDCatLineFluxes_all(verbose=True):
    """
    Run get_LSDCatLineFluxes() for all ids in MiG1D inspections

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    cec.get_LSDCatLineFluxes_all()

    """
    path = '/Users/kschmidt/work/MUSE/ciii_candidates/'
    dat1 = np.genfromtxt(path+'MiG1D_CIIIemitter_candidates_within3DHST_inspection160804_final.txt',
                         dtype=None,names=True,skip_header=2)
    dat2 = np.genfromtxt(path+'MiG1D_CIIIemitter_candidates_withinMUSE_inspection160804_final.txt',
                         dtype=None,names=True,skip_header=2)
    dat3 = np.genfromtxt(path+'MiG1D_CIVemitter_candidates_withinMUSE_inspection160804_final.txt',
                         dtype=None,names=True,skip_header=2)

    MUSEids = dat1['ID'].tolist()+dat2['ID'].tolist()+dat3['ID'].tolist()
    cec.get_LSDCatLineFluxes(MUSEids,outfile='./LSDCatLineFluxes_all_F_3KRON.txt',fluxcol='F_3KRON')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_LSDCatLineFluxes(MUSEids,outfile='./test_F_3KRON.txt',fluxcol='F_3KRON',
                         catparentdir='/Users/kschmidt/work/catalogs/MUSE_GTO/original_per_field_v2.1/',
                         redshiftcat='/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits',
                         verbose=True):
    """
    Pull out the LSDcat line fluxes from the original_per_field catalogs.

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    cec.get_LSDCatLineFluxes([11503085,10414050,10213086],outfile='./test_F_3KRON.txt')

    """
    Nobj = len(MUSEids)
    if verbose: print ' - Will pull out the line fluxes for the '+str(Nobj)+' MUSE IDs provided'

    if verbose: print ' - Getting redshifts from redshift catalog: \n   '+redshiftcat
    dat_z     = pyfits.open(redshiftcat)[1].data
    redshifts = []
    for mi in MUSEids:
        redshifts.append(dat_z['REDSHIFT'][dat_z['UNIQUE_ID'] == str(mi)])

    if verbose: print ' - Will match the found fluxes to the location of: '
    linelist       = MiGs.linelistdic(listversion='rest_uv_main')
    Nlinelistlines = len(linelist)
    if verbose:
        for ll in [linelist[key][0] for key in linelist.keys()]:
            print '   '+str(ll)

    if verbose: print ' - Prepare output ascii file: '+outfile
    fout = open(outfile,'w')
    fout.write('# Linefluxes matched to expected location of lines \n')
    fout.write('# Extracted the fluxes in column '+fluxcol+' of the cataogs from '+catparentdir+' \n')
    fout.write('# Fluxes correspond to detected line fluxes of nearest line. \n')
    fout.write('# lam_match_* gives the distance in wavelength to the location of the flux \n')
    fout.write('# id  redshift ')
    for key in linelist.keys():
        fout.write('     f_'+key+' ferr_'+key+' lam_'+key+' lam_match_'+key)
    fout.write(' \n')

    if verbose: print ' - Looping over object, pulling out fluxes and matching to expected position of lines given z'
    for ii, id in enumerate(MUSEids):
        objz    = redshifts[ii]
        fieldid = str(id)[1:3]
        globres = glob.glob(catparentdir+'cat_ident_candels-cdfs-'+fieldid+'_rid_fluxes*')
        if len(globres) == 0:
            if verbose: print ' - WARNING: did not find any catalog for field '+fieldid+' skipping ID = '+str(id)
        elif len(globres) > 1:
            if verbose: print ' - WARNING: found more than 1 catalog for field '+fieldid+' skipping ID = '+str(id)
        else:
            dat_cat = pyfits.open(globres[0])[1].data
            lineent = np.where(dat_cat['UNIQUE_ID'] == str(id))[0]
            Nlines  = len(lineent)
            if verbose:
                infostr = '   Found '+str(Nlines)+' lines for ID = '+str(id)+' to match to the '+str(Nlinelistlines)+\
                          ' lines in the linelist '
                sys.stdout.write("%s\r" % infostr)
                sys.stdout.flush()

            outputstr = str(id)+' '+str("%.4f" % objz)+' '

            for key in linelist.keys():
                wave_exp      = linelist[key][1]*(1.0+objz) # expected position of line in linelist give obj redshift

                linepositions = dat_cat['LAMBDA_PEAK_FLUX'][lineent]
                dwaves        = linepositions - wave_exp
                bestent       = lineent[ np.where( np.abs(dwaves) == np.min(np.abs(dwaves)) )[0] ]

                line_flux    = dat_cat[fluxcol][bestent]
                line_fluxerr = dat_cat[fluxcol+'_ERR'][bestent]
                line_lam     = dat_cat['LAMBDA_PEAK_FLUX'][bestent]
                line_dlam    = dat_cat['LAMBDA_PEAK_FLUX'][bestent]-wave_exp

                outputstr = outputstr+'    '+\
                            str("%8.2f" % line_flux)+' '+str("%8.2f" % line_fluxerr)+' '+\
                            str("%8.2f" % line_lam )+' '+str("%8.2f" % line_dlam)+'  '

            fout.write(outputstr+'\n')
    if verbose: print '\n - Done...'
    fout.close()
    if verbose: print ' - Wrote output to '+outfile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def print_LSDCatLineFluxSelection(LSDCatFluxFile,matchtol=50,printLyaInfo=False,verbose=True):
    """
    Print objects with flux measurements of lines in the LSDCat cubes within a given match distance
    in wavelength. Fluxes are extracted with get_LSDCatLineFluxes*()

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    LSDCatFluxFile = './LSDCatLineFluxes_all_F_3KRON_160819.txt'
    cec.print_LSDCatLineFluxSelection(LSDCatFluxFile,matchtol=50)

    """
    if verbose: print ' - Loading LSDCat flux data in '+LSDCatFluxFile
    fluxdat = np.genfromtxt(LSDCatFluxFile,names=True,skip_header=4,comments='#',dtype=None)

    if verbose: print ' - Loading MiG1D inspextions '
    MiG1D_path  = '/Users/kschmidt/work/MUSE/ciii_candidates/'
    MiG1D_1     = MiG1D_path+'MiG1D_CIIIemitter_candidates_within3DHST_inspection160804_final.txt'
    MiG1D_1dat  = np.genfromtxt(MiG1D_1,names=True,skip_header=2,comments='#',dtype=None)
    MiG1D_2     = MiG1D_path+'MiG1D_CIIIemitter_candidates_withinMUSE_inspection160804_final.txt'
    MiG1D_2dat  = np.genfromtxt(MiG1D_2,names=True,skip_header=2,comments='#',dtype=None)
    MiG1D_3     = MiG1D_path+'MiG1D_CIVemitter_candidates_withinMUSE_inspection160804_final.txt'
    MiG1D_3dat  = np.genfromtxt(MiG1D_3,names=True,skip_header=2,comments='#',dtype=None)

    Nobj = len(fluxdat)
    if verbose: print '   Found '+str(Nobj)+' objects in file to go through'

    if verbose: print ' - Will look for good matchets (|dlam| < '+str(matchtol)+'A) to the lines: '
    linelist       = MiGs.linelistdic(listversion='rest_uv_main')
    Nlinelistlines = len(linelist)
    if verbose:
        for ll in [linelist[key][0] for key in linelist.keys()]:
            print '   '+str(ll)

    if verbose: print '#  id       line      f        ferr      lam     lam_match     MIG1Dflag1  MIG1Dflag2  MIG1Dflag3  '
    if verbose: print '#                 [1e-20cgs] [1e-20cgs]  [A]        [A]        ciii3dhst   ciiiMUSE    civMUSE  '

    for oo, objid in enumerate(fluxdat['id']):
        objdat = fluxdat[oo]
        for linekey in linelist.keys():
            if np.abs(objdat['lam_match_'+linekey]) < matchtol:
                if (linekey == 'lya') & (printLyaInfo == False):
                    pass
                else:
                    outputstr = str(objid)+'  '+str("%6s" %linekey)+'  '+\
                                str("%8.2f" % objdat['f_'+linekey])+'  '+\
                                str("%8.2f" % objdat['ferr_'+linekey])+'  '+\
                                str("%8.2f" % objdat['lam_'+linekey])+'  '+\
                                str("%8.2f" % objdat['lam_match_'+linekey])+'  '

                    for MiG1D_dat in [MiG1D_1dat,MiG1D_2dat,MiG1D_3dat]:
                        MiGent    = np.where(MiG1D_dat['ID'] == objid)[0]
                        if len(MiGent) == 0:
                            outputstr = outputstr + ' '+str("%10.f" % -99)+' '
                        elif len(MiGent) > 1:
                            sys.exit('More than one match in MiG1D output to '+str(objid))
                        else:
                            for MiGline in MiG1D_dat.dtype.names[1:-5]:
                                if MiGline.split('_')[0].lower() in linekey:
                                    outputstr = outputstr + ' '+str("%10.f" % int(MiG1D_dat[MiGline][MiGent]))+' '

                    if verbose: print outputstr
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

    cec.plot_MUSElya(10414050,3.6609,voffset=300.0,plotSN=False,yrangefull=[-400,1200],showsky=True,wavetype='vac')
    cec.plot_MUSElya(10414050,3.6609,voffset=300.0,plotSN=True,yrangefull=[-3,20],wavetype='vac')

    id       = 11503085 # 12133078
    redshift = 3.7098   # 4.2646
    cec.plot_MUSElya(id,redshift,voffset=-400.0,plotSN=False,yrangefull=[-400,2200],showsky=True,wavetype='vac')
    cec.plot_MUSElya(id,redshift,voffset=-400.0,plotSN=True,yrangefull=[-3,40],wavetype='vac')

    % --------- bad shading of G141 - generate manually ---------
    cec.plot_MUSElya(11015049,4.1451,voffset=300.0,plotSN=False,yrangefull=[-400,1200],showsky=True,wavetype='vac')
    cec.plot_MUSElya(11015049,4.1451,voffset=300.0,plotSN=True,yrangefull=[-3,20],wavetype='vac')

    cec.plot_MUSElya(12350190,4.4955,voffset=300.0,plotSN=False,yrangefull=[-400,1200],showsky=True,wavetype='vac')
    cec.plot_MUSElya(12350190,4.4955,voffset=300.0,plotSN=True,yrangefull=[-3,20],wavetype='vac')

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

    # plot_Lyg      = False; checkwave = 973.*(1+redshift)
    # if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
    #         ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
    #     plot_Lyg      = True

    plot_Lyb      = False; checkwave = 1026.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_Lyb      = True

    plot_Lya      = True

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

    plot_HeII     = False; checkwave = 1640.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_HeII      = True

    plot_CIII     = False; checkwave = 1908.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_CIII     = True

    plot_CIIb     = False; checkwave = 2326.*(1+redshift)
    # if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
    #         ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
    #     plot_CIIb     = True

    plot_MgII     = False; checkwave = 2795.*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_MgII     = True

    plot_OII      = False; checkwave = 3727.5*(1+redshift)
    if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
            ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
        plot_OII      = True

    plot_Hd       = False; checkwave = 4101.74*(1+redshift)
    # if ((checkwave > wavecov_MUSE[0]) & (checkwave < wavecov_MUSE[1])) or  \
    #         ((checkwave > wavecov_3dhst[0]) & (checkwave < wavecov_3dhst[1])):
    #     plot_Hd       = True

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if showsky:
        if verbose: print ' - Getting sky spectra to plot '
        e24 = True
        if e24:
            if verbose: print ' - Assuming objects comes from e24'
            museinfo = pyfits.open('/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits')[1].data
            objent   = np.where(museinfo['UNIQUE_ID'] == str(MUSEid))[0]
        else:
            if verbose: print ' - Assuming objects comes from e36'
            museinfo = pyfits.open('/Users/kschmidt/work/catalogs/MUSE_GTO/merged_catalog_e36_v1.0.fits')[1].data
            objent   = np.where(museinfo['ID'] == int(MUSEid))[0]
        if len(objent) != 1:
            sys.exit(' Found '+str(len(objent))+' matches to '+str(MUSEid)+' in MUSE info fits file... that is weird!')


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
    if plot_Lyb:
        plt.subplot(4, 3, 1) # Lybeta+OIV
        windowcenter = 1030.0
        windowwidth  = 20.0
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
        xrange_Lyb = xrange
        yrange_Lyb = yrange
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 1)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nLy$\\beta$ + OVI doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_Lya:
        plt.subplot(4, 3, 2) # Lyalpha+NV
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
    else:
        plt.subplot(4, 3, 2)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nLy$\\alpha$ + NV doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CII:
        plt.subplot(4, 3, 3) # CII
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
    else:
        plt.subplot(4, 3, 3)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nCII$\\lambda$1336',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_SiIVOIV:
        plt.subplot(4, 3, 4) # SiIV+OIV
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
    else:
        plt.subplot(4, 3, 4)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nSiIV and OIV] doublets',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CIV:
        plt.subplot(4, 3, 5) # CIV
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
    else:
        plt.subplot(4, 3, 5)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nCIV doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_HeII:
        plt.subplot(4, 3, 6) # HeII
        windowcenter = 1653.0
        windowwidth  = 26
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
        xrange_HeII = xrange
        yrange_HeII = yrange

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.subplot(4, 3, 6)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nHeII and OIII] doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_CIII:
        plt.subplot(4, 3, 7) # CIII]
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
    else:
        plt.subplot(4, 3, 7)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nCIII] doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # if plot_CIIb:
    #     plt.subplot(4, 3, XX) # CIIb
    #     windowcenter = 2326.
    #     windowwidth  = 10
    #     if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
    #     xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale
    #
    #     plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
    #     if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)
    #
    #     try:
    #         fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
    #         fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
    #         yrange  = [fluxmin,1.1*fluxmax]
    #     except:
    #         yrange = [0,10]
    #
    #     if fileexist_3dhst:
    #         for f3 in file_3dhst:
    #             f3_dat   = pyfits.open(f3)[1].data
    #             wave     = f3_dat['WAVE_AIR']
    #             flux     = f3_dat['FLUX']-f3_dat['CONTAM']
    #             fluxerr  = f3_dat['FLUXERR']
    #             if plotSN: flux = flux/fluxerr
    #             contam   = f3_dat['CONTAM']
    #
    #             try:
    #                 fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
    #                 fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
    #                 yrange  = [fluxmin,1.1*fluxmax]
    #             except:
    #                 pass
    #
    #             plt.plot(wave, flux, '-',alpha=0.8,color=speccol)
    #
    #             if not plotSN:
    #                 plt.plot(wave, contam, '--',alpha=0.8,color=speccol)
    #
    #                 filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
    #                 filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
    #                 fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
    #                 fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
    #                 plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)
    #
    #     if showsky:
    #         if windowcenter*(redshift+1.0) < 1e4:
    #             sky_w   = sky_wave_muse
    #             sky_f   = sky_flux_muse
    #         else:
    #             sky_w   = sky_wave_nir
    #             sky_f   = sky_flux_nir
    #
    #         skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
    #         skywave = sky_w[skyent]
    #         skylow  = np.zeros(len(skywave))
    #         skyflux = sky_f[skyent]
    #         skyhigh = skyflux
    #         if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])
    #
    #         plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
    #     # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #     cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
    #     # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #     plt.xlim(xrange)
    #     plt.ylim(yrange)
    #     xrange_CIIb = xrange
    #     yrange_CIIb = yrange
    #
    #     plt.xlabel(xlabel)
    #     plt.ylabel(ylabel)
    # else:
    #     plt.subplot(4, 3, XX)
    #
    #     plt.plot(-1,-1)
    #     plt.text(0.5,0.5,'No good coverage of\nCII]',
    #              color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
    #     plt.yticks([])
    #     plt.xticks([])
    #     plt.ylabel('')
    #     plt.xlabel('')
    #     plt.xlim([0,1])
    #     plt.ylim([0,1])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_MgII:
        plt.subplot(4, 3, 8) # MgII
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
    else:
        plt.subplot(4, 3, 8)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\nMgII doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if plot_OII:
        plt.subplot(4, 3, 9) # [OII]
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
    else:
        plt.subplot(4, 3, 9)

        plt.plot(-1,-1)
        plt.text(0.5,0.5,'No good coverage of\n[OII] doublet',
                 color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
        plt.yticks([])
        plt.xticks([])
        plt.ylabel('')
        plt.xlabel('')
        plt.xlim([0,1])
        plt.ylim([0,1])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # if plot_Hd:
    #     plt.subplot(4, 3, XX) # Hd
    #     windowcenter = 4101.74
    #     windowwidth  = 10
    #     if windowcenter*(redshift+1.0) > 10000.: windowwidth  = 70
    #     xrange       = np.asarray([windowcenter-windowwidth,windowcenter+windowwidth])*(redshift+1)/wavescale
    #
    #     plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color=speccol)
    #     if not plotSN: plt.fill_between(wave_MUSE,filllow_MUSE,fillhigh_MUSE,alpha=0.20,color=speccol)
    #
    #     try:
    #         fluxmin = np.min(np.asarray([0, np.min(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
    #         fluxmax = np.max(np.asarray([10, np.max(flux_MUSE[(wave_MUSE > xrange[0]) & (wave_MUSE < xrange[1])]) ]))
    #         yrange  = [fluxmin,1.1*fluxmax]
    #     except:
    #         yrange = [0,10]
    #
    #     if fileexist_3dhst:
    #         for f3 in file_3dhst:
    #             f3_dat   = pyfits.open(f3)[1].data
    #             wave     = f3_dat['WAVE_AIR']
    #             flux     = f3_dat['FLUX']-f3_dat['CONTAM']
    #             fluxerr  = f3_dat['FLUXERR']
    #             if plotSN: flux = flux/fluxerr
    #             contam   = f3_dat['CONTAM']
    #             try:
    #                 fluxmin = np.min(np.asarray([yrange[0], np.min(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
    #                 fluxmax = np.max(np.asarray([yrange[1], np.max(flux[(wave > xrange[0]) & (wave < xrange[1])]) ]))
    #                 yrange  = [fluxmin,1.1*fluxmax]
    #             except:
    #                 pass
    #
    #             plt.plot(wave, flux, '-',alpha=0.8,color=speccol)
    #
    #             if not plotSN:
    #                 plt.plot(wave, contam, '--',alpha=0.8,color=speccol)
    #
    #                 filllow  = flux[(wave > xrange[0]) & (wave < xrange[1])]-fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
    #                 filllow[filllow < np.min(flux)-1e3]  = yrangefull[0]
    #                 fillhigh = flux[(wave > xrange[0]) & (wave < xrange[1])]+fluxerr[(wave > xrange[0]) & (wave < xrange[1])]
    #                 fillhigh[fillhigh > np.max(flux)+1e3]  = yrangefull[1]
    #                 plt.fill_between(wave[(wave > xrange[0]) & (wave < xrange[1])],filllow,fillhigh,alpha=0.20,color=speccol)
    #
    #     if showsky:
    #         if windowcenter*(redshift+1.0) < 1e4:
    #             sky_w   = sky_wave_muse
    #             sky_f   = sky_flux_muse
    #         else:
    #             sky_w   = sky_wave_nir
    #             sky_f   = sky_flux_nir
    #
    #         skyent  = np.where((sky_w > xrange[0]) & (sky_w < xrange[1]))[0]
    #         skywave = sky_w[skyent]
    #         skylow  = np.zeros(len(skywave))
    #         skyflux = sky_f[skyent]
    #         skyhigh = skyflux
    #         if windowcenter*(redshift+1.0) > 1e4: skyhigh = skyhigh / np.max(skyflux) * (yrange[1]-yrange[0])
    #
    #         plt.fill_between(skywave,skylow+yrange[0],skyhigh+yrange[0],alpha=0.3,color='black')
    #     # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #     cec.plot_lines(voffset,llistdic,wavescale,windowwidth,Fsize,col_linemarker,xrange,yrange,redshift,wavetype,LW)
    #     # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #
    #     plt.xlim(xrange)
    #     plt.ylim(yrange)
    #     xrange_Hd = xrange
    #     yrange_Hd = yrange
    #
    #     plt.xlabel(xlabel)
    #     plt.ylabel(ylabel)
    # else:
    #     plt.subplot(4, 3, XX)
    #
    #     plt.plot(-1,-1)
    #     plt.text(0.5,0.5,'No good coverage of\nH\\delta doublet',
    #              color='red',size=Fsize*2,horizontalalignment='center',verticalalignment='center')
    #     plt.yticks([])
    #     plt.xticks([])
    #     plt.ylabel('')
    #     plt.xlabel('')
    #     plt.xlim([0,1])
    #     plt.ylim([0,1])
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
    if plot_Lyb:
        plt.plot(xrange_Lyb,np.zeros(2)+yrange_Lyb[0],'-',color='black',lw=LW)
        plt.plot(xrange_Lyb,np.zeros(2)+yrange_Lyb[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Lyb[0],yrange_Lyb,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Lyb[1],yrange_Lyb,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_Lyb)),yrange_Lyb[1]+0.03*Dyrangefull,'Ly$\\beta$+OVI',
                 color=col_linemarker,size=Fsize,horizontalalignment='center',verticalalignment='bottom')

    if plot_Lya:
        plt.plot(xrange_Lya,np.zeros(2)+yrange_Lya[0],'-',color='black',lw=LW)
        plt.plot(xrange_Lya,np.zeros(2)+yrange_Lya[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Lya[0],yrange_Lya,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_Lya[1],yrange_Lya,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_Lya)),yrange_Lya[1]+0.03*Dyrangefull,'Ly$\\alpha$+NV',
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

    if plot_HeII:
        plt.plot(xrange_HeII,np.zeros(2)+yrange_HeII[0],'-',color='black',lw=LW)
        plt.plot(xrange_HeII,np.zeros(2)+yrange_HeII[1],'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_HeII[0],yrange_HeII,'-',color='black',lw=LW)
        plt.plot(np.zeros(2)+xrange_HeII[1],yrange_HeII,'-',color='black',lw=LW)
        plt.text(np.mean(np.asarray(xrange_HeII)),yrange_HeII[1]+0.03*Dyrangefull,'HeII+OIII]',rotation='vertical',
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

    interactiveplotting = True
    if interactiveplotting:
        plt.plot(wave_MUSE, flux_MUSE, '-',alpha=0.8,color='black',lw=5)
        plt.show()

        array = np.asarray(flux_MUSE.tolist()*500).reshape(500,flux_MUSE.shape[0])
        plt.imshow(array,cmap='gray'); plt.show()

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

            if ll in ['oiv1','oiv2','ovi1','ovi2']:
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
                           yrangefullflux=[-400,1200],yrangefullSN=[-3,30],wavetype='vac',clobber=True,verbose=True):
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

    import ciiiEmitterCandidates as cec
    MUSEidlist = [137030059,137035086]
    zlist      = [0.5236,1.21875]
    cec.plot_MUSElya_forsample(MUSEidlist,zlist,voffsetlist=300.0,datadir='/Users/kschmidt/work/MUSE/spectra1D/Arche170127/spectra/',outputdir='./')


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
        skip = False
        objz = redshiftlist[oo]
        voff = voffsetlist[oo]

        Nfilesexist = len(glob.glob(outputdir+'*'+str(objID)+'*'))
        if verbose:
            idno    = oo+1
            infostr = '   Plotting object '+str(objID)+' at z = '+str(objz)+' indicating voff = '+str(voff)+\
                      '  ('+str(idno)+'/'+str(Nobj)+')                            '

            if (clobber == False) & (Nfilesexist == 2):
                infostr = infostr.replace(')                         ',
                                          ') --> skip (clobber=False)')
                skip = True

            sys.stdout.write("%s\r" % infostr)
            sys.stdout.flush()

        if not skip:
            cec.plot_MUSElya(objID,objz,voffset=voff,plotSN=False,yrangefull=yrangefullflux,outputdir=outputdir,
                             showsky=True,verbose=False,wavetype=wavetype,datadir=datadir)
            cec.plot_MUSElya(objID,objz,voffset=voff,plotSN=True,yrangefull=yrangefullSN,outputdir=outputdir,
                             verbose=False,wavetype=wavetype,datadir=datadir)

    if verbose: print '\n - Done...'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_EWliteraturecomparison(EWcat='/Users/kschmidt/work/catalogs/EWliteraturecatalog.txt',
                                redshiftcolor=True,plotpath='./',plotOuterSymbol=True,verbose=True):
    """
    Plot EW sample from literature using scrips in GitHub/GLASS/GrismReduction/rxj2248_BooneBalestraSource.py

    --- EXAMPLE OF USE ---
    cec.plot_EWliteraturecomparison()

    """
    if verbose: print ' - Load data from '+EWcat
    ewdat = np.genfromtxt(EWcat,dtype=None,skip_header=27,names=True,comments='#')

    xlabel     = 'EW( Ly$\\alpha$ ) / [\\AA]'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xline      = 'lya'
    yline      = 'ciii'
    plotname   = plotpath+'literatureEWs_'+xline+'VS'+yline+'.pdf'
    if redshiftcolor: plotname = plotname.replace('.pdf','_zcol.pdf')
    ylabel     = 'EW( CIII] ) / [\\AA]'

    goodent    = np.where((ewdat['ew_'+xline] != -9999) & (ewdat['ew_'+yline] != -9999) &
                          (ewdat['ew_'+xline] != -99) & (ewdat['ew_'+yline] != -99))[0]
    if len(goodent) > 0:
        xrange     = [-45,220]#1000]
        yrange     = [-5,80]
        xrangezoom = [-20,40]
        yrangezoom = [-4,12]
        zcut       = [2.8,3.9]

        bbs.plot_EWcatalog(xline,yline,ewdat[goodent],plotname,xlabel,ylabel,
                           xrange,yrange,xrangezoom,yrangezoom,zcut=zcut,
                           plotOuterSymbol=plotOuterSymbol,redshiftcolor=redshiftcolor,verbose=verbose)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xline      = 'lya'
    yline      = 'civ'
    plotname   = plotpath+'literatureEWs_'+xline+'VS'+yline+'.pdf'
    if redshiftcolor: plotname = plotname.replace('.pdf','_zcol.pdf')
    ylabel     = 'EW( CIV ) / [\\AA]'

    goodent    = np.where((ewdat['ew_'+xline] != -9999) & (ewdat['ew_'+yline] != -9999) &
                          (ewdat['ew_'+xline] != -99) & (ewdat['ew_'+yline] != -99))[0]
    if len(goodent) > 0:
        xrange     = [-45,220]#1000]
        yrange     = [-5,45]
        xrangezoom = [-20,60]
        yrangezoom = [-4,-1.5]
        zcut       = [2.8,3.9]

        bbs.plot_EWcatalog(xline,yline,ewdat[goodent],plotname,xlabel,ylabel,
                           xrange,yrange,xrangezoom,yrangezoom,zcut=zcut,
                           plotOuterSymbol=plotOuterSymbol,redshiftcolor=redshiftcolor,verbose=verbose)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def MiG1Doutput2LaTeXtable(MiG1Doutput,match3dhst,verbose=True,
                           MUSEWideInfo='/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits'):
    """
    Convert a MiG1D inspection output file to a latex table (adding extrac information on the sources)
    and print the results to the screen.

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    MiG1Doutput = './MiG1D_CIIIemitter_candidates_withinMUSE_inspection160804_final.txt'
    match3dhst  = './CIIIemitter_candidates_withinMUSE.txt'
    cec.MiG1Doutput2LaTeXtable(MiG1Doutput,match3dhst=match3dhst,verbose=True)

    MiG1Doutput = './MiG1D_CIIIemitter_candidates_within3DHST_inspection160804_final.txt'
    match3dhst  = './CIIIemitter_candidates_within3DHST.txt'
    cec.MiG1Doutput2LaTeXtable(MiG1Doutput,match3dhst=match3dhst,verbose=True)

    MiG1Doutput = './MiG1D_CIVemitter_candidates_withinMUSE_inspection160804_final.txt'
    match3dhst  = './CIVemitter_candidates_withinMUSE.txt'
    cec.MiG1Doutput2LaTeXtable(MiG1Doutput,match3dhst=match3dhst,verbose=True)

    """
    if verbose: print ' - Loading data'
    dat_MiG1D = np.genfromtxt(MiG1Doutput,dtype=None,names=True,comments='#',skip_header=2)
    dat_match = np.genfromtxt(match3dhst ,dtype=None,names=True,comments='#',skip_header=4)
    dat_MWide = pyfits.open(MUSEWideInfo)[1].data

    objIDs = dat_MiG1D['ID']
    Nobj   = len(objIDs)
    if verbose: print ' - Collect data for '+str(Nobj)+' objects in MiG1D output and write table format to screen\n'

    if verbose: print 'ID         & R.A. & Dec.  & ID     & $r_\\textrm{Match}$ & $z$ & ' \
                      'Ly$\\beta$ & OVI & NV & SiIV,OIV & CIV & HeII & CIII & MgII & ' \
                      'Close Neighbor & $z_\\textrm{Offset}$ & $Qz_\\textrm{Offset}$ & v$_\\textrm{Offset}$ \\\\ '

    if verbose: print 'MUSE-Wide  & [Deg]& [Deg] & 3D-HST & [arcsec] & MUSE-Wide & ' \
                      '1026\\AA & 1032,1035\\AA & 1239,1242\\AA & 1294--1403\\AA & 1548,1551\\AA & 1640\\AA & 1907,1909\\AA & 2796,2803\\AA &' \
                      ' & & & [km/s] \\\\ \n\\hline'

    for objid in objIDs:
        latexstr = str(objid)
        ent_MiG1D = np.where(dat_MiG1D['ID'] == objid)[0]
        ent_match = np.where(dat_match['ID_MUSE'] == objid)[0]
        ent_MWide = np.where(dat_MWide['UNIQUE_ID'] == str(objid))[0]

        latexstr  = latexstr+' & '+str("%.8f" % dat_MWide['RA'][ent_MWide])+' & '+str("%.8f" % dat_MWide['DEC'][ent_MWide])
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if len(ent_match) == 1:
            rmatch = dat_match['R_MATCH'][ent_match]
            if rmatch < 0.5:
                latexstr  = latexstr+' & '+str("%.f" % dat_match['ID_3DHST'][ent_match])+' & '+str("%.4f" % rmatch)
            else:
                latexstr  = latexstr+' &  &  '
        else:
            latexstr  = latexstr+' &  &  '
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        latexstr  = latexstr+' & '+str("%.4f" % dat_MWide['REDSHIFT'][ent_MWide])+'$\\pm$'+\
                    str("%.4f" % dat_MWide['REDSHIFT_ERR'][ent_MWide])
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        MiG1Dstr =  ' & '+str("%.f" % dat_MiG1D['Lyb_1025'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['OVI_1032'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['NV_1240'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['SiIVOIV_1399'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['CIV_1549'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['HeII_1640'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['CIII_1909'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['MgII_2800'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['close_counterpart'][ent_MiG1D])+\
                    ' & '+str("%.4f" % dat_MiG1D['byhandredshift'][ent_MiG1D])+\
                    ' & '+str("%.f" % dat_MiG1D['byhandredshift_quality'][ent_MiG1D])+' '

        MiG1Dstr  = MiG1Dstr.replace('& 0 ','& - ').replace('& -99.0000 ','&  ')
        latexstr  = latexstr+MiG1Dstr
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        z_bh     = dat_MiG1D['byhandredshift'][ent_MiG1D]
        if z_bh != -99.0 :
            z_lya    = dat_MWide['REDSHIFT'][ent_MWide]
            cc       = 299792.458 # km/s
            voffset  = (z_lya - z_bh) * cc / (z_bh + 1.0)
            latexstr  = latexstr+' & '+str("%.4f" % voffset)
        else:
            latexstr  = latexstr+' & '
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose: print latexstr+' \\\\'
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def measurelinefluxes(MUSEids,outputdir='./',generatelinelists=True,measurefluxes=True,lines_manual={},
                      MUSEcat='/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits',
                      clobber=False,verbose=True,verbose_flux=False,SNthreshinit=1.0, SNthreshanal=1.0,
                      skipLSDCatIfExists = False):
    """

    Wrapper for import fluxmeasurementsMUSEcubes functions to measure line fluxes and limits in MUSE
    cubes at expected locations of rest-frame UV lines incl CIII and CIV

    --- INPUT ---
    MUSEids              List of object IDs (UNIQUE_ID) to estimate lines fluxes for
    outputdir            Directory to save generated line lists and flux files to
    generatelinelists    Set to True to generate line list. If False the line list are assumed to
                         exist in outputdirectory
    measurefluxes        Set to True to measure line fluxes. If False the flux outputs are assumed to
                         exist in outputdirectory
    lines_manual         Provide a dictionary with entries for 'lines' (cube positions) to estimate
                         flux at. Expects the format:
                         lines_manual = {MUSEID:[ [LINENAME1,LINENAME2,...], [WAVE_OBS1, WAVE_OBS2,...] ]}
    MUSEcat              Muse catalog to get redshift,ra and dec of object in MUSEids list from
    clobber              If True outputs will be overwritten
    SNthreshinit         Set S/N threshold for LSDCat flux determination
    SNthreshanal         Set S/N analysis threshold for LSDCat flux determination
    skipLSDCatIfExists   To enable conitnuing a full run skip the timeconsuming LSDCat step if ouput already exists.
    verbose              Toggle verbosity of main script
    verbose_flux         Toggle verbosity of flux functions

    --- EXAMPLE OF USE ---

    import ciiiEmitterCandidates as cec
    MUSEids  = ['11503085','10306046']
    manlines = {'11503085':[ ['testman1','testman2'], [7501.99, 7502.99] ]}
    fluxcats = cec.measurelinefluxes(MUSEids,lines_manual = manlines)

    """
    fielddic = {'1':'cdfs','2':'cosmos'}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading obj info from MUSE catalog:\n   '+MUSEcat
    MUSEdat = pyfits.open(MUSEcat)[1].data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print " - Grabbing list of emission lines from MiGs' linelist"
    linesall  = MiGs.linelistdic(listversion='full')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nobj = str(len(MUSEids))
    if verbose: print ' - Looping over the '+Nobj+' MUSE IDs provided '
    fluxcatalogs = []
    for oo, obj_id in enumerate(MUSEids):
        if verbose: print ' ------------------ Object '+obj_id+' ('+str(oo+1)+'/'+Nobj+') ------------------'
        if verbose: print ' - Resetting list of line locations to MiG lines'
        linewaves = []
        linenames = []
        for line in linesall.keys():
            linenames.append(linesall[line][0])
            linewaves.append(linesall[line][1])
        linenames = np.asarray(linenames)
        linewaves = np.asarray(linewaves)

        MUSEent  = np.where(MUSEdat['UNIQUE_ID'] == obj_id)[0]
        if len(MUSEent) != 1:
            sys.exit(' Something is not right - found '+str(len(MUSEent))+' entries matching MUSEid='+obj_id+' in the catalog '+MUSEcat)

        linecat  = outputdir+'/'+obj_id+'_linelist.fits'
        obj_z    = MUSEdat['REDSHIFT'][MUSEent][0]
        obj_ra   = MUSEdat['RA'][MUSEent][0]
        obj_dec  = MUSEdat['DEC'][MUSEent][0]

        wave_obs = linewaves*(obj_z+1) # NB: Halpha outside MUSE range

        if obj_id in lines_manual.keys():
            if verbose: print ' - Found line positions in lines_manual; adding them to default line list'
            man_names = lines_manual[obj_id][0]
            man_waves = lines_manual[obj_id][1]
            for mm in xrange(len(man_names)):
                linenames    = np.append(linenames,man_names[mm])
                man_wave_obs = man_waves[mm]
                wave_obs     = np.append(wave_obs,man_wave_obs)

        objIDs   = np.asarray( [obj_id]*len(wave_obs) )
        lineIDs  = np.asarray( [obj_id+str("%.3d" % (nn+1)) for nn in xrange(len(wave_obs))] )
        ras      = np.asarray( [obj_ra]*len(wave_obs) )
        decs     = np.asarray( [obj_dec]*len(wave_obs) )

        pointing = obj_id[1:3]
        fieldno  = obj_id[0]
        cube     = '/Volumes/DATABCKUP2/MUSE-Wide/data/'+\
                   '/median_filtered_DATACUBE_candels-'+fielddic[fieldno]+'-'+pointing+'_v1.0.fits_effnoised_dcbgc.fits'

        if os.path.isfile(cube):
            if generatelinelists:
                if verbose: print ' - Generating line list'
                if (skipLSDCatIfExists == True) & (os.path.isfile(linecat)):
                    if verbose: print '   ... but skipLSDCatIfExists=True and '+linecat+' exists so moving on'
                else:
                    fmm.save_LSDCatFriendlyFitsFile(linecat,lineIDs,objIDs,ras,decs,wave_obs,radecwave=True,coordcube=cube,
                                                    linenames=linenames,clobber=clobber,verbose=verbose_flux)
            else:
                if verbose: print ' - Skip generating line list - assume it exists'

            if measurefluxes:
                nowstr  = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
                if verbose: print ' - Measuring fluxes for line list (takes ~3 min per object; now it is '+nowstr+') '
                output_expect = linecat.replace('.fits','_fluxes.fits')
                if (skipLSDCatIfExists == True) & (os.path.isfile(output_expect)):
                    if verbose: print '   ... but skipLSDCatIfExists=True and '+output_expect+' exists so moving on'
                    fluxcatalog = output_expect
                else:
                    fluxcatalog = fmm.measure_fluxes(linecat, field=fielddic[fieldno], field_id=pointing,
                                                     verbose=verbose_flux, clobber=clobber,
                                                     SNthreshinit=SNthreshinit, SNthreshanal=SNthreshanal)
                    nowstr  = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
                    if verbose: print '   ... done at '+nowstr
            else:
                if verbose: print ' - Skipping measuring line fluxes - assume output exists'
                fluxcatalog = linecat.replace('.fits','_fluxes.fits')

            fluxcatalogs.append(fluxcatalog)
        else:
            if verbose: print '   WARNING: Did not find the filtered cube:\n   '+cube
            fluxcatalogs.append(outputdir+'/'+obj_id+'_NoCubeFoundForObject')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    return fluxcatalogs
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def measurelinefluxes_allobj(C3inMUSE=True,C3inbetween=False,C3in3DHST=False,measurefluxes=True,
                             generatelinelists=True,verbose=True,verbose_flux=False,clobber=False,
                             SNthreshinit=1.0, SNthreshanal=1.0,skipLSDCatIfExists=False):
    """
    Wrapper for measurelinefluxes() listing the IDs to measure lines fluxes for and defining dictionary of any
    manual line positions.

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    fluxcats = cec.measurelinefluxes_allobj(C3inMUSE=True,C3inbetween=False,C3in3DHST=False,measurefluxes=True)

    fluxcats = cec.measurelinefluxes_allobj(C3inMUSE=True,C3inbetween=True,C3in3DHST=True,measurefluxes=True)

    """
    MUSEids_combined = cec.getids(C3inMUSE=C3inMUSE,C3inbetween=C3inbetween,C3in3DHST=C3in3DHST)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    lines_manual = {}
    lines_manual ['11503085'] = [ ['testman1','testman2'], [7501.99, 7502.99] ]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fluxcats = cec.measurelinefluxes(MUSEids_combined,lines_manual = lines_manual,clobber=clobber,
                                     measurefluxes=measurefluxes,generatelinelists=generatelinelists,
                                     verbose=verbose,verbose_flux=verbose_flux,SNthreshinit=SNthreshinit,
                                     SNthreshanal=SNthreshanal,skipLSDCatIfExists=skipLSDCatIfExists)

    return fluxcats

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def getids(C3inMUSE=True,C3inbetween=False,C3in3DHST=False):

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # CIII in MUSE
    MUSEids_combined = []
    if C3inMUSE:
        MUSEids = ['10105016','10109024','10110025','10111026','10112027','10212085','10213086','10214087','10215088',
                   '10219097','10301001','10302002','10305045','10306046','10308055','10309056','10310057','10311058',
                   '10404029','10407038','10409041','10412048','10414050','10415052','10505027','10507031','10509034',
                   '10603018','10604019','10607033','10608034','10609035','10612043','10614046','10615047','10704077',
                   '10708086','10715098','10716099','10806104','10809107','10810108','10903027','10904028','10906034',
                   '10908036','10909037','10910038','11002003','11003005','11004006','11006008','11007009','11101001',
                   '11102002','11205037','11206038','11207039','11208040','11305024','11307030','11406067','11503085',
                   '11504089','11508100','11512104','11520114','11521115','11522116','11524120','11528125','11602059',
                   '11605062','11611072','11615086','11616087','11619092','11621094','11701001','11702002','11705019',
                   '11706020','11802011','11806033','11901001','11902002','11903003','11904004','11905005','11906006',
                   '11907007','12003003','12004004','12006008','12101001','12103013','12104014','12107020','12112031',
                   '12113032','12202034','12204055','12205056','12207063','12210066','12211067','12212068','12213069',
                   '12221111','12222112','12306090','12307091','12316117','12318120','12331149','12401001','12403012',
                   '12404013','12405014']
        #MUSEids = ['12133078','11503085','10414050','10806104','10509034']
        MUSEids_combined = MUSEids_combined + MUSEids
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # CIII between MUSE and 3D-HST
    if C3inbetween:
        MUSEids = ['10115031','10116032','10123043','10124044','10125045','10318081','10319082','10326097','10331105',
                   '10332106','10333107','10419059','10424069','10619055','10621059','10625070','10627072','10632083',
                   '10724120','10817130','10914056','10915057','10917059','11015049','11018062','11109022','11110023',
                   '11113028','11115030','11116031','11214060','11314051','11410076','11532151','11535156','11536157',
                   '11623096','11624097','11711043','11717058','11814057','11818066','11821073','11916033','11917034',
                   '11923056','11926065','11927066','12009015','12011018','12117043','12118044','12119045','12133078',
                   '12228123','12232127','12235132','12348186','12350190','12351191','12353193','12411026','12414030',
                   '12415031','12424054','12426056','12427057','12428058']
        MUSEids_combined = MUSEids_combined + MUSEids
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # CIII in 3D-HST
    if C3in3DHST:
        MUSEids = ['10127047','10130050','10432086','10635088','10733150','10734151','10741159','10830156','10924083',
                   '11121036','11222074','11419101','11639142','11728077','11729078','11739091','11834094','11929068',
                   '11931070','11936075','11938077','11940079','12124068','12132077','12135080','12236133','12361203']
        MUSEids_combined = MUSEids_combined + MUSEids
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return MUSEids_combined
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_fluxmeasurements(objects,fluxcatdir='./',outname='./fluxcomparison.pdf',showids=False,verbose=True,
                          line2plot='Lya',fluxcolumn='F_KRON',xrange=[-10000,25000],yrange=[-10000,25000],
                          LSDCat='/Users/kschmidt/work/catalogs/MUSE_GTO/candels_1-24_emline_master_v2.1.fits',
                          tdhstMatch='/Users/kschmidt/work/MUSE/candelsCDFS_3DHST/MUSECDFS_z0p0-7p0_cmtol10p0_v2p1_3DHSTinfo.fits',tdhstcol=False,tdhstMatchTol=0.5,markobj=False):
    """
    Compare measured fluxes (and other quantities) for a list of objects.

    --- INPUT ---
    objects
    fluxcatdir
    outname
    showids
    verbose
    line2plot       Name of line to compare/look for in flux catalogs, i.e. from the LSDCat
                    original_per_field_v2.1 outputs
    fluxcolumn      Column name of flux columns to compare
    LSDCat
    tdhstMatch      Catalog containing match of MUSE-Wide sources to the 3D-HST data
    tdhstcol        If a column name from the 3D-HST catalogs is provided, this will be plotted on the x-axis.
                    Use a list of [clon_name,error_name] to specify what to plot.
    tdhstMatchTol   Tolerance on matches to 3D-HST catalogs.
    markobj         Provide list of objects to mark their points with red color for easy recognition

    --- EXAMPLE OF USE ---
    import ciiiEmitterCandidates as cec
    MUSEids = cec.getids(C3inMUSE=True,C3inbetween=False,C3in3DHST=False)
    cec.plot_fluxmeasurements(MUSEids,outname='./fluxcomp_Lya_1KRON_C3inMUSEids.pdf',fluxcolumn='F_KRON',fluxcatdir='./ForceFluxC3inMUSE_fullrun161031/',xrange=[-500,15000],yrange=[-500,15000])
    cec.plot_fluxmeasurements(MUSEids,outname='./fluxcomp_Lya_2KRON_C3inMUSEids.pdf',fluxcolumn='F_2KRON',fluxcatdir='./ForceFluxC3inMUSE_fullrun161031/',xrange=[0,25000],yrange=[-5000,25000])
    cec.plot_fluxmeasurements(MUSEids,outname='./fluxcomp_Lya_3KRON_C3inMUSEids.pdf',fluxcolumn='F_3KRON',fluxcatdir='./ForceFluxC3inMUSE_fullrun161031/',xrange=[0,25000],yrange=[-5000,25000],showids=True)


    cec.plot_fluxmeasurements(MUSEids,outname='./fluxcomp_CIV_1KRON_C3inMUSEids.pdf',line2plot='CIV',fluxcolumn='F_KRON',fluxcatdir='./ForceFluxC3inMUSE_fullrun161031/',xrange=[-500,15000],yrange=[-500,15000])


    cec.plot_fluxmeasurements(MUSEids,outname='./fluxcomp_Lya_1KRON_C3inMUSEids_Mass.pdf',line2plot='Lya',fluxcolumn='F_KRON',fluxcatdir='./ForceFluxC3inMUSE_fullrun161031/',tdhstcol=['lmass'],xrange=[7,12],yrange=[-500,15000],showids=True)



    """
    fielddic = {'1':'cdfs','2':'cosmos'}
    if verbose: print ' - Loading data to compare objects with:'
    if verbose: print '   LSDCat overview    : '+LSDCat
    lsdcatdat  = pyfits.open(LSDCat)[1].data
    if verbose: print '   3D-HST match info : '+tdhstMatch
    tdhstdat = pyfits.open(tdhstMatch)[1].data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plotname = outname
    if showids: plotname = plotname.replace('.pdf','_wIDs.pdf')
    if verbose: print ' - Setting up plot '+plotname
    fig = plt.figure(figsize=(6, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.15, right=0.95, bottom=0.10, top=0.95)
    Fsize    = 10
    lthick   = 1
    marksize = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    if tdhstcol:
        plt.title('Match between MUSE-Wide and 3D-HST $<$ '+str(tdhstMatchTol)+'arcsec',fontsize=Fsize)

    plt.plot([np.min(xrange+yrange),np.max(xrange+yrange)],
             [np.min(xrange+yrange),np.max(xrange+yrange)],
             ls='--',color='k',lw=lthick)


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    linenamedic = {}
    linenamedic['Lya']    = 'Ly$\\alpha$ $\\lambda$1216'
    # linenamedic['lyb']    = 'Ly$\\beta$ $\\lambda$1025'
    # linenamedic['lyg']    = 'Ly$\gamma$ $\\lambda$973'
    # linenamedic['nv1']    = 'NV $\\lambda$1239'
    # linenamedic['nv2']    = 'NV $\\lambda$1243'
    # linenamedic['cii']    = 'CII $\\lambda$1336'
    # linenamedic['Siiv1']  = 'SiIV $\\lambda$1394'
    # linenamedic['Siiv2']  = 'SiIV $\\lambda$1403'
    # linenamedic['oiv1']   = 'OIV $\\lambda$1397'
    # linenamedic['oiv2']   = 'OIV $\\lambda$1400'
    linenamedic['CIV']   = 'CIV $\\lambda$1548'
    # linenamedic['civ2']   = 'CIV $\\lambda$1551'
    # linenamedic['heii']   = 'HeII $\\lambda$1640'
    # linenamedic['oiiib1'] = 'OIII] $\\lambda$1661'
    # linenamedic['oiiib2'] = 'OIII] $\\lambda$1666'
    # linenamedic['ciii1']  = '[CIII] $\\lambda$1907'
    # linenamedic['ciii2']  = 'CIII] $\\lambda$1909'
    # linenamedic['ciib']   = 'CII] $\\lambda$2326'
    # linenamedic['mgii1']  = 'MgII] $\\lambda$2796'
    # linenamedic['mgii2']  = 'MgII] $\\lambda$2803'
    # linenamedic['oii1']   = '[OII] $\\lambda$3726'
    # linenamedic['oii2']   = '[OII] $\\lambda$3729'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for objid in objects:
        if tdhstcol:
            tdhstent    = np.where(tdhstdat['id_MUSE'] == objid)[0]
            tdhstid     = tdhstdat['id_3DHST'][tdhstent]
            tdhstRmatch = tdhstdat['r_match_arcsec'][tdhstent]
            if tdhstRmatch <= tdhstMatchTol:
                plottdhst = True
            else:
                if verbose: print ' - Poor match (tdhstMatchTol='+str(tdhstMatchTol)+') to 3D-HST catalog for '+objid
                plottdhst = False

        fluxcat  = False
        fcsearch = fluxcatdir+objid+'_linelist_fluxes.fits'
        if os.path.isfile(fcsearch): fluxcat = fcsearch
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        pointcol = 'black'
        if markobj:
            if objid in markobj:
                pointcol = 'red'

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if fluxcat:
            forcedat    = pyfits.open(fcsearch)[1].data
            forceent    = np.where(forcedat['LINENAME'] == linenamedic[line2plot])[0]

            if len(forceent) == 1:
                forceline    = forcedat['LINENAME'][forceent][0]
                forceflux    = forcedat[fluxcolumn][forceent]
                forcefluxerr = forcedat[fluxcolumn+'_ERR'][forceent]

                fieldno     = objid[0]
                fieldid     = objid[1:3]
                fieldname   = 'candels-'+fielddic[fieldno]+'-'+str(fieldid)
                LSDCatpath  = '/Users/kschmidt/work/catalogs/MUSE_GTO/original_per_field_v2.1/'
                LSDfieldcat = LSDCatpath+'cat_ident_'+fieldname+'_rid_fluxes_v2.1.fits'
                LSDfluxdat  = pyfits.open(LSDfieldcat)[1].data
                LSDCobjent  = np.where((LSDfluxdat['UNIQUE_ID'] == objid) & (LSDfluxdat['Identification'] == line2plot))[0]

                if tdhstcol:
                    if plottdhst:
                        tdhstval = tdhstdat[tdhstcol[0]][tdhstent]
                        try:
                            tdhsterr = tdhstdat[tdhstcol[1]][tdhstent]
                        except:
                            tdhsterr = None

                        if verbose: print ' - Plotting quanteties for object '+objid
                        if np.isfinite(tdhstval):
                            dxval = np.abs(xrange[1]-xrange[0])*0.05
                            dyval = np.abs(yrange[1]-yrange[0])*0.05
                            if tdhstval < xrange[0]:
                                plt.arrow(xrange[0]+dxval, forceflux[0], -dxval*0.66, 0,
                                          head_width=dyval*0.33, head_length=dxval*0.33, fc=pointcol, ec=pointcol)

                                xvalplot = xrange[0]+dxval
                                haval    = 'left'
                            elif tdhstval > xrange[1]:
                                plt.arrow(xrange[1]-dxval, forceflux[0], dxval*0.66, 0,
                                          head_width=dyval*0.33, head_length=dxval*0.33, fc=pointcol, ec=pointcol)

                                xvalplot = xrange[1]-dxval
                                haval    = 'right'
                            else:
                                plt.errorbar(tdhstval,forceflux,xerr=tdhsterr,yerr=forcefluxerr,color=pointcol,
                                             ls='o',lw=lthick)
                                xvalplot = tdhstval
                                haval    = 'left'

                            if showids:
                                plt.text(xvalplot,forceflux,str(objid),color=pointcol,ha=haval,va='bottom')
                else:
                    if len(LSDCobjent) == 1:
                        LSDCatflux    = LSDfluxdat[fluxcolumn][LSDCobjent]
                        LSDCatfluxerr = LSDfluxdat[fluxcolumn+'_ERR'][LSDCobjent]

                        if verbose: print ' - Plotting fluxes for object '+objid
                        plt.errorbar(LSDCatflux,forceflux,xerr=LSDCatfluxerr,yerr=forcefluxerr,color=pointcol,
                                     ls='o',lw=lthick)

                        if showids:
                            plt.text(LSDCatflux,forceflux,str(objid),color=pointcol,ha='left',va='bottom')
                    else:
                        if verbose: print '   No entry in LSDCat (looked for '+objid+' and '+line2plot+') --> Skipping'
            else:
                if verbose: print '   No entry in force flux catalog (looked for '+linenamedic[line2plot]+') --> Skipping'
        else:
            if verbose: print '   No force flux catalog (looked for '+fcsearch+') --> Skipping'


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plt.xlim(xrange)
    plt.ylim(yrange)

    if tdhstcol:
        tdhstnamedic = {}
        tdhstnamedic['lmass']  = r'log10(M$_{\rm{FAST 3D-HST}}$ / [Msun])'
        tdhstnamedic['lsfr']   = r'log10(SFR$_{\rm{FAST 3D-HST}}$ / [Msun/yr])'
        tdhstnamedic['lage']   = r'log10(Age$_{\rm{FAST 3D-HST}}$ / [yr])'
        #tdhstnamedic['sfr']   = 'log10(SFR$_\rm{W14 3D-HST}$/[Msun/yr])'
        tdhstnamedic['jh_mag'] = r'M$_{\rm{F140W 3D-HST}}$ / [AB Mag]'

        plt.xlabel(tdhstnamedic[tdhstcol[0]])

        if tdhstcol[0] in ['jh_mag']:
            plt.xlim(xrange[::-1])
    else:
        plt.xlabel(r'LSDCat flux; '+fluxcolumn.replace('_','\_')+'('+line2plot+') [1e-20 erg/s/cm$^2$]')

    plt.ylabel(r'Force flux; '+fluxcolumn.replace('_','\_')+'('+forceline+') [1e-20 erg/s/cm$^2$]')

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


    # # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # if verbose: print ' - Plotting EWs for object '+objid
    # if EWcat:
    #     if verbose: print '   EW plotting not enabled yet --> Skipping'
    # else:
    #     if verbose: print '   No flux catalog found (looked for '+ecsearch+') --> Skipping'
    #
    #
    # # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def calculateEWs4objlist(objects,fluxcatdir='./',verbose=True):
    """
    Compare measured fluxes (and other quanteties) for a list of objects.

    --- EXAMPLE OF USE ---

    """

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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =