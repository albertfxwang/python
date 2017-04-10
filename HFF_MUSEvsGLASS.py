# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import numpy as np
import pyfits
import sys
import pdb
import crossmatch as cm
import matplotlib.pyplot as plt
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def crossmatch_catalogs(verbose=True):
    """
    Command to crossmatch catalogs...

    """
    # ------------------------------
    # -----------A2744--------------
    # ------------------------------
    catin      = '/Users/kschmidt/work/GLASS/CRALvisit161012/A2744_conf1.fits'
    colin      = ['#IDsec','RA','DEC']
    namein ='MUSEa2744conf1'

    catmatch   = '/Users/kschmidt/work/catalogs/GLASScatalogs/GLASScatalog_A2744_150515.fits'
    colmatch   = ['NUMBER','X_WORLD','Y_WORLD']
    namematch='GLASSv001maincat'
    # ------------------------------
    catin      = '/Users/kschmidt/work/GLASS/CRALvisit161012/A2744_conf1.fits'
    colin      = ['#IDsec','RA','DEC']
    namein='MUSEa2744conf1'

    catmatch   = '/Users/kschmidt//work/GitHub/GLASS/docs/redshiftcatalogs/hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_redshiftcatalog.txt'
    colmatch   = ['0','1','2']
    namematch='GLASSv001zcat'
    # ------------------------------
    # -----------M0416--------------
    # ------------------------------
    catin      = '/Users/kschmidt/work/GLASS/CRALvisit161012/MACS0416_z_highconfidence.fits'
    colin      = ['ID','RA','DEC']
    namein ='MUSEm0416highconf'

    catmatch   = '/Users/kschmidt/work/catalogs/GLASScatalogs/GLASScatalog_MACS0416.1-2403_150515.fits'
    colmatch   = ['NUMBER','X_WORLD','Y_WORLD']
    namematch='GLASSv001maincat'
    # ------------------------------
    catin      = '/Users/kschmidt/work/GLASS/CRALvisit161012/MACS0416_z_highconfidence.fits'
    colin      = ['ID','RA','DEC']
    namein ='MUSEm0416highconf'

    catmatch   = '/Users/kschmidt//work/GitHub/GLASS/docs/redshiftcatalogs/hlsp_glass_hst_wfc3_macs0416-fullfov-pa999_ir_v001_redshiftcatalog.txt'
    colmatch   = ['0','2','3']
    namematch='GLASSv001zcat'
    # ------------------------------
    # -----------COMMAND------------
    # ------------------------------
    outputfile = './matchCat_'+namein+'_'+namematch+'.txt'
    matcharr   = cm.matchCat(catin,colin,catmatch,colmatch,namein=namein,namematch=namematch,outputfile=outputfile)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def plot_redshifts(cmcat,zcat_muse,zcat_glass,r_match_tol=0.5,Nhdrcols_glass=30,
                   idcol_muse='ID',idcol_glass='ID',zcol_muse='z',zcol_glass='redshift',Qzcol_glass='redshift_quality',
                   plotname='./redshiftplot.pdf',zrange=[-1.1,6.7],dzmark=0.25,show_IDs=False,verbose=True):
    """

    Plot the redshift with matches in both MUSE and GLASS catalogs
    --- INPUT ---
    verbose         Toggle verbosity

    --- EXAMPLE OF USE ---
    import HFF_MUSEvsGLASS as hmg

    cmcat      = './matchCat_MUSEa2744conf1_GLASSv001zcat.fits'
    zcat_muse  = '/Users/kschmidt/work/GLASS/CRALvisit161012/A2744_conf1.fits'
    zcat_glass = '/Users/kschmidt//work/GitHub/GLASS/docs/redshiftcatalogs/hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_redshiftcatalog.txt'
    plotname   = './redshiftcomparison_MUSEa2744conf1_GLASSv001zcat.pdf'

    hmg.plot_redshifts(cmcat,zcat_muse,zcat_glass,plotname=plotname,idcol_muse='#IDsec',Nhdrcols_glass=30,show_IDs=True)

    cmcat      = './matchCat_MUSEa2744full_GLASSv001zcat.fits'
    zcat_muse  = '/Users/kschmidt/work/GLASS/CRALvisit161012/A2744_redshifts_cat_final.fits'
    zcat_glass = '/Users/kschmidt//work/GitHub/GLASS/docs/redshiftcatalogs/hlsp_glass_hst_wfc3_a2744-fullfov-pa999_ir_v001_redshiftcatalog.txt'
    plotname   = './redshiftcomparison_MUSEa2744full_GLASSv001zcat.pdf'

    hmg.plot_redshifts(cmcat,zcat_muse,zcat_glass,plotname=plotname,idcol_muse='#ID',zcol_muse='Z',Nhdrcols_glass=30,show_IDs=True)




    cmcat      = './matchCat_MUSEm0416highconf_GLASSv001zcat.fits'
    zcat_muse  = '/Users/kschmidt/work/GLASS/CRALvisit161012/MACS0416_z_highconfidence.fits'
    zcat_glass = '/Users/kschmidt//work/GitHub/GLASS/docs/redshiftcatalogs/hlsp_glass_hst_wfc3_macs0416-fullfov-pa999_ir_v001_redshiftcatalog.txt'
    plotname   = './redshiftcomparison_MUSEm0416highconf_GLASSv001zcat.pdf'

    hmg.plot_redshifts(cmcat,zcat_muse,zcat_glass,plotname=plotname,idcol_muse='ID',idcol_glass='id',Nhdrcols_glass=16,show_IDs=True,zcol_muse='Redshift',zcol_glass='zgrism')

    cmcat      = 'matchCat_MUSEm0416SWhighconf_GLASSv001zcat.fits'
    zcat_muse  = '/Users/kschmidt/work/GLASS/CRALvisit161012/MACS0416-SW-high-condidence.fits'
    zcat_glass = '/Users/kschmidt//work/GitHub/GLASS/docs/redshiftcatalogs/hlsp_glass_hst_wfc3_macs0416-fullfov-pa999_ir_v001_redshiftcatalog.txt'
    plotname   = './redshiftcomparison_MUSEm0416SWhighconf_GLASSv001zcat.pdf'

    hmg.plot_redshifts(cmcat,zcat_muse,zcat_glass,plotname=plotname,idcol_muse='ID',idcol_glass='id',Nhdrcols_glass=16,show_IDs=True,zcol_muse='Z',zcol_glass='zgrism')


    """
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Loading data in MUSE and 3DHST catalogs'
    dat_MUSE  = pyfits.open(zcat_muse)[1].data
    dat_GLASS = np.genfromtxt(zcat_glass,names=True,dtype=None,skip_header=Nhdrcols_glass,comments='#')

    id_MUSEcat     = dat_MUSE[idcol_muse].astype(float)
    id_GLASScat    = dat_GLASS[idcol_glass].astype(float)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Get redshift for matchest within the match tolerance '+str(r_match_tol)+' arcsec'
    dat_cmatch = pyfits.open(cmcat)[1].data

    goodmatch  = np.where(dat_cmatch[dat_cmatch.dtype.names[6]] <= r_match_tol)[0]
    Ngood      = len(goodmatch)

    if Ngood == 0:
        sys.exit('No crosmatches within tolerance of '+str(r_match_tol)+' arcsec')
    else:
        if verbose: print '   Will plot redshift for the '+str(Ngood)+' objects satisfying the match tolerance '
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose: print ' - Setting up and generating plot'
    fig = plt.figure(figsize=(6, 5))
    fig.subplots_adjust(wspace=0.1, hspace=0.1,left=0.1, right=0.95, bottom=0.10, top=0.95)
    Fsize    = 14
    lthick   = 1
    marksize = 7
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=Fsize)
    plt.rc('xtick', labelsize=Fsize)
    plt.rc('ytick', labelsize=Fsize)
    plt.clf()
    plt.ioff()
    plt.title('d$z$\_mark= '+str(dzmark)+'; r\_match\_tol='+str(r_match_tol)+' arc sec',fontsize=Fsize)

    plt.plot(zrange,zrange,'--k',lw=lthick)

    if show_IDs: print ' id_MUSE  z_MUSE  id_GLASS  z_GLASS  dz  Qz_GLASS'


    for objent in goodmatch:
        id_MUSE  = dat_cmatch[dat_cmatch.dtype.names[0]][objent]
        id_GLASS = dat_cmatch[dat_cmatch.dtype.names[3]][objent]

        MUSEent  = np.where(id_MUSEcat == id_MUSE)[0]
        zMUSE    = float(dat_MUSE[zcol_muse][MUSEent][0])

        GLASSent = np.where(id_GLASScat == id_GLASS)[0]
        zGLASS   = float(dat_GLASS[zcol_glass][GLASSent][0])
        QzGLASS  = float(dat_GLASS[Qzcol_glass][GLASSent][0])

        dz = np.abs(zMUSE - zGLASS)

        if dz < dzmark:
            color = 'green'
        else:
            color = 'red'

        if zGLASS == -1:
            color = 'blue'

        plt.errorbar(zMUSE,zGLASS,xerr=None,yerr=None,fmt='o',lw=lthick,
                     ecolor=color, markersize=marksize,markerfacecolor=color,markeredgecolor = 'k')

        if show_IDs:
            infostring = str("%10s" % id_MUSE)+' '+str("%8.4f" % zMUSE)+' '+str("%10s" % id_GLASS)+' '+\
                         str("%8.4f" % zGLASS)+' '+str("%8.5f" % dz)+str("%6.2f" % QzGLASS)

            if zGLASS >= 0:
                print ' '+infostring
            else:
                pass
                #print '#'+infostring

            #plt.text(zMUSE,zGLASS,' idMUSE='+str(id_MUSE)+', idGLASS='+str(id_GLASS),fontsize=Fsize)

    plt.xlabel(r'$z$ MUSE', fontsize=Fsize)
    plt.ylabel(r'$z$ GLASS', fontsize=Fsize)

    plt.xlim(zrange)
    plt.ylim(zrange)

    #--------- LEGEND ---------
    # plt.errorbar(-1,-1,xerr=None,yerr=None,fmt='o',lw=lthick,ecolor='blue', markersize=marksize*2,
    #              markerfacecolor='blue',markeredgecolor = 'k',label='EA$z$Y Photo')
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