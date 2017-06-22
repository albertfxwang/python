# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import pdb
import os
import commands
import sys
import glob
import pyfits
import numpy as np
import fits2ascii as f2a
import kbsutilities as kbs
import crossmatch as cm
import MUSE_AOsetup as mao
from astropy import units
import collections
from astropy.vo.client import conesearch
from astropy.coordinates import ICRS, FK5
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def clusterinfo(verbose=True):
    """
    Return dictionary with useful cluster information

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    cldic = mao.clusterinfo()
    """

    cluster_dic = {}
    cluster_dic['A370']            = {'ra':39.97,'dec': -1.57667,'name':'A370'}
    cluster_dic['A2744']           = {'ra':3.59750,'dec':-30.39056,'name':'A2744'}
    cluster_dic['MACS0416.1-2403'] = {'ra':64.03913,'dec': -24.06775,'name':'MACS0416'}
    cluster_dic['MACS0717.5+3745'] = {'ra':109.38187,'dec': 37.75514,'name':'MACS0717'}
    cluster_dic['MACS0744.9+3927'] = {'ra':116.22,'dec': 39.45678,'name':'MACS0744'}
    cluster_dic['MACS1149.6+2223'] = {'ra':177.39942,'dec': 22.39861,'name':'MACS1149'}
    cluster_dic['MACS1423.8+2404'] = {'ra':215.949,'dec': 24.07792,'name':'MACS1423'}
    cluster_dic['MACS2129.4-0741'] = {'ra':322.35858,'dec': -7.69133,'name':'MACS2129'}
    cluster_dic['RXJ1347.5-1145']  = {'ra':206.87746,'dec': -11.75281,'name':'RXJ1347'}
    cluster_dic['RXJ2248']         = {'ra':342.18458,'dec': -44.52667,'name':'RXJ2248'}

    return cluster_dic
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def MUSE_AOregion_DeepCall(filename,ra,dec,color='magenta',verbose=True):
    """
    Generating region file for AO setup for Deep Field MUSE call
    http://www.eso.org/sci/activities/docs/Call_MUSE_DF.pdf

    --- INPUT ---

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    cluster_dic = mao.clusterinfo()
    for key in cluster_dic.keys():
        cli     = cluster_dic[key]
        region  = mao.MUSE_AOregion_DeepCall('MUSE_AO_DS9region_'+cli['name']+'.reg',cli['ra'],cli['dec'])

    """
    green_box_dist = 1.0/3600.0*37.0 # 37 arcsec offset
    green_box_N    = [ra,dec+green_box_dist]
    green_box_S    = [ra,dec-green_box_dist]
    green_box_E    = [ra+green_box_dist,dec]
    green_box_W    = [ra-green_box_dist,dec]

    regionstring = """
# Region file format: DS9 version 4.1
fk5
# composite(%s,%s,0) || composite=1 color=%s width=3
box(%s,%s,50",50",0) || # color=yellow width=2
box(%s,%s,30",8",0) || # color=green width=2
box(%s,%s,30",8",0) || # color=green width=2
#box(%s,%s,8",30",0) || # color=green width=2
#box(%s,%s,8",30",0) || # color=green width=2
box(%s,%s,60",60",0) || # color=red width=2
circle(%s,%s,52") || # color=blue width=4
circle(%s,%s,101.5") || # color=blue width=4
circle(%s,%s,0.5") # color=%s width=4
"""  % (ra,dec,color,
        ra,dec,
        green_box_N[0],green_box_N[1],
        green_box_S[0],green_box_S[1],
        green_box_E[0],green_box_E[1],
        green_box_W[0],green_box_W[1],
        ra,dec,
        ra,dec,
        ra,dec,
        ra,dec,color)

    fout = open(filename,'w')
    fout.write(regionstring)
    fout.close()
    if verbose: print ' - Generated region file '+filename
    return regionstring
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def perform_conesearch(ra,dec,search_radius,catalog=None,verbose=True):
    """
    Perform a conesearch using astropy's VO capabilities

    --- INPUT ---
    ra                  richt ascension [deg]
    dec                 declination     [deg]
    search_radius       radius around ra and dec to search [deg]
    catalo              name of catalog to search; set to None to get option printed
    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    catalog         = 'guide*star'
    catalog         = 'Guide Star Catalog v2 1'
    search_result   = mao.perform_conesearch(177.3992,22.39885,0.05,catalog=catalog)
    """

    if catalog == None:
        print ' - No catalog selected, please choose from the following list (wildcards ok, e.g., "guide*star":'
        for cat in conesearch.list_catalogs(): print '  ',cat
        return None
    else:
        if verbose: print ' - Performing cone search in radius ',search_radius,\
            ' deg around (ra,dec) = (',ra,',',dec,')'

    cat_list = conesearch.list_catalogs(pattern=catalog)

    obj_coord = ICRS(ra=ra*units.deg, dec=dec*units.deg)

    search_result = conesearch.search_all(obj_coord, search_radius*units.deg, catalog_db=cat_list)

    if verbose:
        print ' - Searched ',len(search_result),' catalogs. Each catalog contains the following number of entries:'
        print ' '
        for url, tab in search_result.items():
            print '   CatUrl: ',url
            print '   TabCol: ',tab.nrows
            print '   Nrows : ',tab.array.data.dtype.names
            print ' '

    return search_result
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def save_conesearch(conesearch_result,filename='Default',emptystring='NoneKBS',verbose=True):
    """
    Save output from perform_conesearch to ascii file

    --- INPUT ---
    conesearch_results   output from perfomr_conesearch
    filename_ext         extension appended ascii filenames

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    mao.save_conesearch(search_result,filename='./testfilename')

    """

    for url, tab in conesearch_result.items():
        if filename == 'Default':
            outputfile   = './' + (str(url.split('/')[-1])).replace('&','amp').replace('.','p')+'.txt'
        else:
            outputfile = filename+'.txt'
        if verbose: print '\n - Saving data table to '+outputfile
        hdr          = 'Cone search results. \nNrows = '+str(tab.nrows)+'. \nurl = '+url
        colnames     = tab.array.data.dtype.names
        hdr          = hdr+'\n'+'  '.join(colnames)
        data         = tab.array.data.copy()

        if verbose: print ' - Checking output array for empty entries (adding "'+emptystring+'" if empty)'
        for col in colnames:
            for vv, val in enumerate(data[col]):
                if val == '':
                    data[col][vv] = emptystring

        # outformat = []
        # for name in colnames:
        #     fmttype = tab.array.data.dtype[name]
        #     if fmttype == float:
        #         outformat.append("%16.8f")
        #     if fmttype == int:
        #         outformat.append("%16i")
        #     if fmttype == 'O':
        #         outformat.append("%16i")
        #     if fmttype == 'S10': #  What about strings in general?
        #         outformat.append("%.50s")

        np.savetxt(outputfile,data,header=hdr,delimiter='  ',fmt='%s')
        if verbose: print ' - Successfully saved file'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def search_GSPCv2p4(outputfile,ra,dec,search_radius,clobber=False,verbose=True):
    """
    Searching Guide Star Photometric Catalog V2p4.fits

    --- INPUT ---
    outputname
    ra
    dec
    search_radius
    clobber
    verbose

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    mao.search_GSPCv2p4('MACS1149test.fits',177.3992,22.39885,0.5,clobber=False,verbose=True)

    """
    catalog = '/Users/kschmidt/work/catalogs/GuideStarPhotometricCatalogV2p4.fits'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    GSPC_hdu      = pyfits.open(catalog) # Load the FITS hdulist
    GSPC_dat      = GSPC_hdu['II_272_gspc24'].data
    GSPC_hdr      = GSPC_hdu['II_272_gspc24'].header

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    GSPC_IDall    = GSPC_dat['ID']
    GSPC_RAall    = GSPC_dat['RAJ2000']
    GSPC_Decall   = GSPC_dat['DEJ2000']
    Nobj_GSPC     = len(GSPC_IDall)
    if verbose: print ' - Will find matches to the '+str(Nobj_GSPC)+' in \n   '+catalog
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    r_match  = np.sqrt( (np.cos(np.rad2deg(dec))*(GSPC_RAall-ra))**2.0 + (GSPC_Decall-dec)**2.0 )

    goodmatch = np.where(r_match <= search_radius)[0]
    if len(goodmatch) > 0:
        if verbose: print ' - Writing output to '+outputfile
        if (os.path.isfile(outputfile)) & (clobber==False): # check if file already exists
            if verbose: print '   WARNING: output file already exists. clobber=False so did not overwrite '
        else:
            if (os.path.isfile(outputfile)) & (clobber==True) & verbose:
                print '   Output file already exists but clobber=True so overwriting it'

            pyfits.writeto(outputfile, GSPC_dat[goodmatch], clobber=clobber)
    else:
        if verbose: print ' - WARNING No matches found within search_radius = '+str(search_radius)+\
                          ' No file returned'
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def search_CLASHcat(cluster,output='default',search_radius=0.08333,clobber=False,verbose=True):
    """

    --- INPUT ---
    cluster        (full) Name of cluster to search
    output         Name of output to generate
    search_radius  Radius to return objects in around cluster ra and dec (from clusterinfo())
    clobber        Overwrite outputfile if it already exists?
    verbose        Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    cluster   = 'MACS2129.4-0741'
    outputcat = mao.search_CLASHcat(cluster,verbose=True)


    """
    cluster_short = cluster.split('.')[0].lower()
    catalog = '/Users/kschmidt/work/catalogs/CLASHcatalogs/hlsp_clash_hst_ir_'+cluster_short+'_cat.txt.FITS'
    ra, dec = mao.clusterinfo()[cluster.upper()]['ra'], mao.clusterinfo()[cluster.upper()]['dec']
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CLASH_hdu      = pyfits.open(catalog) # Load the FITS hdulist
    CLASH_dat      = CLASH_hdu[1].data
    CLASH_hdr      = CLASH_hdu[1].header

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CLASH_IDall    = CLASH_dat['id']
    CLASH_RAall    = CLASH_dat['ra']
    CLASH_Decall   = CLASH_dat['dec']
    Nobj_CLASH     = len(CLASH_IDall)
    if verbose: print ' - Will find matches to the '+str(Nobj_CLASH)+' in \n   '+catalog
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if output== 'default':
        outputfile = cluster_short.upper()+'_CLASHcatSearchOutput_rsearch'+\
                     str(search_radius).replace('.','p')+'.fits'
    else:
        outputfile = output
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    r_match  = np.sqrt( (np.cos(np.deg2rad(dec))*(CLASH_RAall-ra))**2.0 + (CLASH_Decall-dec)**2.0 )

    goodmatch = np.where(r_match <= search_radius)[0]
    if len(goodmatch) > 0:
        if verbose: print ' - Writing output to '+outputfile
        if (os.path.isfile(outputfile)) & (clobber==False): # check if file already exists
            if verbose: print '   WARNING: output file already exists. clobber=False so did not overwrite '
        else:
            if (os.path.isfile(outputfile)) & (clobber==True) & verbose:
                print '   Output file already exists but clobber=True so overwriting it'

            pyfits.writeto(outputfile, CLASH_dat[goodmatch], clobber=clobber)
    else:
        if verbose: print ' - WARNING No matches found within search_radius = '+str(search_radius)+\
                          ' No file returned'

    return outputfile

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def search_ROMANcat(cluster,output='default',search_radius=0.08333,clobber=False,verbose=True):
    """

    --- INPUT ---
    cluster        (full) Name of cluster to search
    output         Name of output to generate
    search_radius  Radius to return objects in around cluster ra and dec (from clusterinfo())
    clobber        Overwrite outputfile if it already exists?
    verbose        Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    cluster   = 'a2744'
    outputcat = mao.search_ROMANcat(cluster,verbose=True)

    """
    cluster_short = cluster.split('.')[0].lower()

    if cluster_short == 'a2744':
        catfile = '/Users/kschmidt/work/GitHub/GLASS/ROMAN_CATALOGS/A2744/A2744_CLUSTER.cat'
    elif cluster_short == 'macs0416':
        catfile = '/Users/kschmidt/work/GitHub/GLASS/ROMAN_CATALOGS/M0416/M0416_CLUSTER.cat'
    elif cluster_short == 'macs1149':
        #coordfile = '/Users/kschmidt/work/GitHub/GLASS/ROMAN_CATALOGS/M1149/outcord.1007_HST.dat'
        #catfile   = '/Users/kschmidt/work/GitHub/GLASS/ROMAN_CATALOGS/M1149/outmag.1007_HST.dat'
        catfile    = '/Users/kschmidt/work/catalogs/ROMANphotocat_M1149combcat.cat'
    else:
        print ' WARNING No ROMAN catalog exists for the cluster '+cluster+'. Nothing returned.'
        return None
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    outpath     = '/Users/kschmidt/work/catalogs/'
    fitscatalog = f2a.ascii2fits(catfile,asciinames=True,skip_header=0,outpath=outpath,verbose=verbose)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ROMAN_hdu      = pyfits.open(fitscatalog) # Load the FITS hdulist
    ROMAN_dat      = ROMAN_hdu[1].data
    ROMAN_IDall    = ROMAN_dat['ID']
    ROMAN_RAall    = ROMAN_dat['RA']
    ROMAN_Decall   = ROMAN_dat['DEC']
    Nobj_ROMAN     = len(ROMAN_IDall)
    if verbose: print ' - Will find matches to the '+str(Nobj_ROMAN)+' in \n   '+catfile
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ra, dec = mao.clusterinfo()[cluster.upper()]['ra'], mao.clusterinfo()[cluster.upper()]['dec']
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if output== 'default':
        outputfile = cluster_short.upper()+'_ROMANcatSearchOutput_rsearch'+\
                     str(search_radius).replace('.','p')+'.fits'
    else:
        outputfile = output
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    r_match  = np.sqrt( (np.cos(np.deg2rad(dec))*(ROMAN_RAall-ra))**2.0 + (ROMAN_Decall-dec)**2.0 )

    goodmatch = np.where(r_match <= search_radius)[0]
    if len(goodmatch) > 0:
        if verbose: print ' - Writing output to '+outputfile
        if (os.path.isfile(outputfile)) & (clobber==False): # check if file already exists
            if verbose: print '   WARNING: output file already exists. clobber=False so did not overwrite '
        else:
            if (os.path.isfile(outputfile)) & (clobber==True) & verbose:
                print '   Output file already exists but clobber=True so overwriting it'

            pyfits.writeto(outputfile, ROMAN_dat[goodmatch], clobber=clobber)
    else:
        if verbose: print ' - WARNING No matches found within search_radius = '+str(search_radius)+\
                          ' No file returned'

    return outputfile
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def create_objectDS9reg(filename,ra_name,dec_name,mag_names,mag_limits=[17.5,20.0,20.5],
                        id_name='objID',skip_header=3,star_name=None,mag_min=5,verbose=True):
    """
    Generate DS9 region files for star (object) catalog with color coding according to magnitude

    --- INPUT ---
    filename         ascii file with object information. E.g., output from save_conesearch
                     Can also take fits binary table input. Checks for '.fit' in filename
    ra_name          Name of righ ascension column
    dec_name         Name of declination column
    mag_names        Names of column with manitudes corresponding to the limits in mag_limits
    mag_limits       Magnitude limits to use
    id_name          Name of ID column
    skip_header      Number of header lines to skip if using np.genfromtxt to load ascii file
    star_name        Name of column containing information about object stellarity
                     Note that the selection is based on the column name, i.e., 'stel'
                     assumes a CLASH-like setup.
    mag_min          The brightest magnitude to consider. Avoids spurious magnitude entries like -99
    verbose          toggle vebosity

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    mao.create_objectDS9reg('testfilename_test160113.txt','ra','dec',['RMag','RMag','VMag'])

    mao.create_objectDS9reg('USNOB1p0_objects.txt','ra','dec',['Rmag','Rmag','Vmag'],skip_header=5,id_name='ID',star_name='Gal')

    mao.create_objectDS9reg('MACS1149test.fits','RAJ2000','DEJ2000',['Rmag','Rmag','Vmag'],id_name='ID')

    mao.create_objectDS9reg('MACS2129_CLASHcatSearchOutput_rsearch0p08333.fits','ra','dec',['f606w_mag','f606w_mag','f606w_mag'],id_name='id',star_name='stel')

    """
    mag_limits   = np.array(mag_limits)
    if type(mag_names) == str: mag_names = [mag_names]*len(mag_limits)
    mag_colors   = ['cyan','yellow','green'] # color ref: http://www.eso.org/sci/activities/docs/Call_MUSE_DF.pdf
    mag_size     = [0.5,0.6,0.7]
    if verbose: print ' - Will generate DS9 region file for objects brighter than ',mag_limits
    if verbose: print '   using color coding ',mag_colors
    if verbose: print '   and sizes          ',mag_size


    if verbose: print ' - Loading data and iterating over objects '
    if '.fits' in filename:
        hdu  = pyfits.open(filename)
        data = hdu[1].data
        reg_filename = filename.replace('.fits','_brightobj.reg')
    else:
        data     = np.genfromtxt(filename,skip_header=skip_header,names=True,comments='#')
        reg_filename = filename.replace('.txt','_brightobj.reg')

    fout = open(reg_filename,'w')

    reg_hdr = """# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
"""

    fout.write(reg_hdr)
    for objdat in data:
        obj_id  = objdat[id_name]
        obj_ra  = objdat[ra_name]
        obj_dec = objdat[dec_name]

        for ii, mag_lim in enumerate(mag_limits):
            mag_ent = np.where(np.asarray(data.dtype.names) == mag_names[ii])
            if np.size(mag_ent) == 0:
                sys.exit(' - Did not find the magnitude column "'+mag_names[ii]+'" --> ABORTING')
            else:
                mag_ent = mag_ent[0][0]

            ObjIsStar = True
            if star_name != None: # checking if star flag indicates a galaxy or a star
                if star_name == 'stel':
                    obj_star = objdat[star_name]
                    if obj_star < 0.8: ObjIsStar = False
                elif star_name == 'Gal':
                    obj_ent  = np.where(data[id_name] == obj_id)[0]
                    data_str = np.genfromtxt(filename,skip_header=skip_header,names=True,comments='#',
                                             dtype='f,40a,f,f,f,f,f,f,f,f,f,f,40a')
                    obj_star = data_str[obj_ent][star_name]
                    if obj_star != 'N': ObjIsStar = False
                else:
                    print 'WARNING No setup for star_name (stellarity flag) "'+star_name+'";',
                    print ' Assuming all objects to be stars '

            obj_mag  = objdat[mag_ent]
            obj_text = str(objdat[id_name])+' '+mag_names[ii]+' = '+str("%.2f" % obj_mag)
            if (obj_mag < mag_lim) and (obj_mag > mag_min):
                reg_col   = mag_colors[ii]
                reg_size  = mag_size[ii]

                obj_row   = 'circle('+str(obj_ra)+','+str(obj_dec)+','+str(reg_size)+'") # color='+\
                            str(reg_col)+' width=2 '
                if ObjIsStar:
                    obj_row = obj_row+' font="times 16 bold roman" text={'+obj_text+'}\n'
                else:
                    obj_row = obj_row+' font="times 16 normal italic" text={'+obj_text+'} dash=1\n'

                fout.write(obj_row)
    fout.close()
    if verbose: print ' - Wrote region file of brigh objects to '+reg_filename

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def open_fits_and_regions(cluster='all',printcommand=False,regstring='',verbose=True):
    """
    Open full FOV files in DS9 and overlay the regions to inspect

    --- INPUT ---
    cluster        Name of cluster to display. Can use 'all' to open all 10 clusters
    printcommand   Set to True to only print the command so you can open it manually from another terminal
    regstring      Only open region files with given string in name (string after the cluster name)
    verbose        Toggle verbosity

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    mao.open_fits_and_regions('MACS1149')
    mao.open_fits_and_regions('A2744',regstring='EDIT',printcommand=True)
    mao.open_fits_and_regions('RXJ2248',regstring='EDIT',printcommand=True)
    mao.open_fits_and_regions('MACS0416',regstring='EDIT',printcommand=True)
    mao.open_fits_and_regions('MACS1149',regstring='EDIT',printcommand=True)

    """

    imgpath = '/Users/kschmidt/work/GLASS/images_fullfov/'
    regpath = '/Users/kschmidt/work/MUSE/call_PublicDeepFied/'

    ds9cmd = ' ds9 '

    if cluster == 'A370' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*A370*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_A370-IR_drz.fits -region '+regfiles

    if cluster == 'A2744' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*A2744*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_hlsp_frontier_hst_wfc3-60mas_abell2744_f140w_v1.0_drz.fits -region '+regfiles

    if cluster == 'MACS0416' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*MACS0416*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_hlsp_frontier_hst_wfc3-60mas_macs0416_f140w_v1.0_drz.fits -region '+regfiles

    if cluster == 'MACS0717' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*MACS0717*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_MACS0717-IR_drz.fits -region '+regfiles

    if cluster == 'MACS0744' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*MACS0744*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_hlsp_clash_hst_wfc3ir_macs0744_total_v1_drz.fits -region '+regfiles

    if cluster == 'MACS1149' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*MACS1149*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_hlsp_clash_hst_wfc3ir_macs1149_total_v1_drz.fits -region '+regfiles

    if cluster == 'MACS1423' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*MACS1423*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_hlsp_clash_hst_wfc3ir_macs1423_total_v1_drz.fits -region '+regfiles

    if cluster == 'MACS2129' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*MACS2129*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_hlsp_clash_hst_wfc3ir_macs2129_total_v1_drz.fits -region '+regfiles

    if cluster == 'RXJ1347' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*RXJ1347*'+regstring+'*.reg'))
        ds9cmd = ds9cmd+' '+imgpath+'refimage_hlsp_clash_hst_wfc3ir_rxj1347_total_v1_drz.fits -region '+regfiles

    if cluster == 'RXJ2248' or cluster == 'all':
        regfiles = ' -region '.join(glob.glob(regpath+'*RXJ2248*'+regstring+'*.reg'))
        ds9cmd   = ds9cmd+' '+imgpath+'refimage_hlsp_clash_hst_wfc3ir_rxj2248_total_v1_drz.fits -region '+regfiles

    ds9cmd = ds9cmd+' '

    if printcommand:
        print ds9cmd+' & '
    else:
        cmdout = commands.getoutput(ds9cmd)

        if cmdout != '':
            print 'Error with command: ds9cmd'
            print ' commands.getoutput output:\n',cmdout
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def run_all_GLASSclusters(search_radius=0.16666,printcommand=True,catalog=None,verbose=True):
    """
    Perform cone search and generate the object region files for all 10 GLASS clusters

    --- INPUT ---
    search_radius       radius around ra and dec to search [deg]
    printcommand        Print DS9 command instead of opening the files?

    --- EXAMPLE OF USE ---
    import MUSE_AOsetup as mao
    mao.run_all_GLASSclusters(0.16666)

    """
    search_catalog     = 'Guide Star Catalog v2 1'
    cluster_dic        = mao.clusterinfo()
    if verbose: print ' - Search catalog set ('+search_catalog+') and cluster info loaded '

    if verbose: print ' - Searching for objects within '+str(search_radius)+' degrees of the cluster '
    for key in cluster_dic.keys():
        cli            = cluster_dic[key]
        if verbose: print '========================='+cli['name']+'========================='
        # -------------------- SEARCH via VO --------------------
        search_result  = mao.perform_conesearch(cli['ra'],cli['dec'],search_radius,
                                                catalog=search_catalog,verbose=verbose)
        conesearchfile = '/Users/kschmidt/work/MUSE/call_PublicDeepFied/'+\
                         cli['name']+'_conesearch_'+search_catalog.replace(' ','_')+\
                         '_Rsearch_'+str(search_radius).replace('.','p')
        mao.save_conesearch(search_result,filename=conesearchfile,verbose=verbose)
        if os.path.isfile(conesearchfile+'.txt'):
            mao.create_objectDS9reg(conesearchfile+'.txt','ra','dec',['RMag','RMag','VMag'],verbose=verbose)

        # -------------------- SEARCH GSPC --------------------
        GSPC_search_radius = 1.0 #deg
        GSPCsearchfile = '/Users/kschmidt/work/MUSE/call_PublicDeepFied/'+\
                         cli['name']+'_GSPCv2p4search_Rsearch_'+str(GSPC_search_radius).replace('.','p')+'.fits'

        mao.search_GSPCv2p4(GSPCsearchfile,cli['ra'],cli['dec'],GSPC_search_radius,clobber=False,verbose=True)
        if os.path.isfile(GSPCsearchfile):
            mao.create_objectDS9reg(GSPCsearchfile,'RAJ2000','DEJ2000',['Rmag','Rmag','Vmag'],id_name='ID')

        # -------------------- SEARCH CLASH --------------------
        if ('A370' not in key) and ('A2744' not in key):
            outputcat = mao.search_CLASHcat(key,search_radius=search_radius,verbose=True)
            mao.create_objectDS9reg(outputcat,'ra','dec',['f606w_mag','f606w_mag','f606w_mag'],
                                    id_name='id',star_name='stel')

        # -------------------- SEARCH ROMAN --------------------
        outputcat = mao.search_ROMANcat(key,search_radius=search_radius,verbose=True)
        if outputcat != None:
            mao.create_objectDS9reg(outputcat,'RA','DEC',['V606','V606','V606'],id_name='ID')

    if printcommand:
        if verbose: print '\n\n - To open the fits files use the DS9 command:'
    mao.open_fits_and_regions('all',printcommand=printcommand,verbose=verbose)
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_extinctions4GLASSclusters(filter='F125W',valreturn='median',radius=1.7/60.,
                                  stepsize=[0.4/60.,0.4/60.],clusters='GLASS',verbose=True):
    """

    Using kbs.getAv_area to get extinction of GLASS clusters

    mao.get_extinctions4GLASSclusters()

    """
    if clusters == 'GLASS':
        cldic = mao.clusterinfo()
    elif clusters == 'others':
        cldic = {}
        cldic['bullet'] = {'ra':104.62958,'dec': -55.94694 , 'name':'Bullet'}
        cldic['A1689']  = {'ra':197.87292,'dec': -1.33806  , 'name':'A1689'}

    for key in cldic.keys():
        if verbose: print ' ---------------- Get extinction for '+cldic[key]['name']+' ----------------'
        if verbose == 'Full':
            vb_getAv_area = True
        else:
            vb_getAv_area = False
        A, EBV, grid = kbs.getAv_area(cldic[key]['ra'],cldic[key]['dec'],radius,stepsize=stepsize,
                                      valreturn=valreturn,filter=filter,verbose=vb_getAv_area)
        if verbose: print ' - The results is A, EBV, grid = ',A,',', EBV,',', grid,','
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def get_extinctions4superGLASSclusters(filter='F125W',valreturn='median',radius=1.7/60.,
                                       stepsize=[0.6/60.,0.6/60.],clusters='GLASS',verbose=True):
    """

    Using kbs.getAv_area to get extinction of GLASS clusters

    Avalarr = mao.get_extinctions4superGLASSclusters()

    """
    cluster_dic = collections.OrderedDict()
    # cluster_dic['RXJ1532      ']            = {'ra':233.22408,'dec':   30.34964  ,'name':'RXJ1532      '}
    # cluster_dic['M0429        ']            = {'ra': 67.40042,'dec':   -2.88556  ,'name':'M0429        '}
    # cluster_dic['A1423        ']            = {'ra':179.32192,'dec':   33.61039  ,'name':'A1423        '}
    # cluster_dic['M1311        ']            = {'ra':197.75696,'dec':   -3.17764  ,'name':'M1311        '}
    # cluster_dic['A209         ']            = {'ra': 22.96904,'dec':  -13.61078  ,'name':'A209         '}
    # cluster_dic['MS2137       ']            = {'ra':325.06325,'dec':  -23.66131  ,'name':'MS2137       '}
    # cluster_dic['A611         ']            = {'ra':120.23679,'dec':   36.05669  ,'name':'A611         '}
    # cluster_dic['CL1226       ']            = {'ra':186.74321,'dec':   33.54650  ,'name':'CL1226       '}
    # cluster_dic['M0329        ']            = {'ra': 52.42367,'dec':   -2.19658  ,'name':'M0329        '}
    # cluster_dic['RXJ2129      ']            = {'ra':322.41642,'dec':    0.08856  ,'name':'RXJ2129      '}
    # cluster_dic['M1720        ']            = {'ra':260.07063,'dec':   35.60656  ,'name':'M1720        '}
    # cluster_dic['M1931        ']            = {'ra':292.95692,'dec':  -26.57611  ,'name':'M1931        '}
    # cluster_dic['A383         ']            = {'ra': 42.01400,'dec':   -3.52908  ,'name':'A383         '}
    # cluster_dic['M0647        ']            = {'ra':101.95846,'dec':   70.24714  ,'name':'M0647        '}
    # cluster_dic['M1206        ']            = {'ra':181.55117,'dec':   -8.80067  ,'name':'M1206        '}
    # cluster_dic['A2261        ']            = {'ra':260.61354,'dec':   32.13294  ,'name':'A2261        '}
    # cluster_dic['M1115        ']            = {'ra':168.96688,'dec':    1.49906  ,'name':'M1115        '}
    # cluster_dic['Bullet       ']            = {'ra':104.65792,'dec':  -55.95000  ,'name':'Bullet       '}
    # cluster_dic['A1689        ']            = {'ra':197.89250,'dec':   -1.36556  ,'name':'A1689        '}
    # cluster_dic['cl1232.5-1250']            = {'ra':188.12708,'dec':  -12.84333  ,'name':'cl1232.5-1250'}
    # cluster_dic['cl1216.8-1201']            = {'ra':184.18792,'dec':  -12.02167  ,'name':'cl1216.8-1201'}
    # cluster_dic['cl1138.2-1133']            = {'ra':174.54292,'dec':  -11.56056  ,'name':'cl1138.2-1133'}
    # cluster_dic['cl1354.2-1230']            = {'ra':208.54042,'dec':  -12.51694  ,'name':'cl1354.2-1230'}
    # cluster_dic['cl1054.4-1146']            = {'ra':163.60208,'dec':  -11.77222  ,'name':'cl1054.4-1146'}
    # cluster_dic['cl1227.9-1138']            = {'ra':186.99542,'dec':  -11.58694  ,'name':'cl1227.9-1138'}
    # cluster_dic['cl1054.7-1245']            = {'ra':163.68167,'dec':  -12.76444  ,'name':'cl1054.7-1245'}
    # cluster_dic['cl1040.7-1155']            = {'ra':160.16833,'dec':  -11.93444  ,'name':'cl1040.7-1155'}
    cluster_dic['SDSSJ0004-0103']           = {'ra':1.216458	,'dec':-1.054389  ,'name':'SDSSJ0004-0103'}
    cluster_dic['SDSSJ0108+0623']           = {'ra':17.175000	,'dec':6.412222   ,'name':'SDSSJ0108+0623'}
    cluster_dic['SDSSJ0150+2725']           = {'ra':27.504167	,'dec':27.426667  ,'name':'SDSSJ0150+2725'}
    cluster_dic['SDSSJ0146-0929']           = {'ra':26.733333	,'dec':-9.497778  ,'name':'SDSSJ0146-0929'}
    cluster_dic['SDSSJ0333-0651']           = {'ra':53.270833	,'dec':-6.856111  ,'name':'SDSSJ0333-0651'}
    cluster_dic['SDSSJ0851+3331']           = {'ra':132.912500	,'dec':33.519722  ,'name':'SDSSJ0851+3331'}
    cluster_dic['SDSSJ1110+6459']           = {'ra':167.575000	,'dec':64.996389  ,'name':'SDSSJ1110+6459'}
    cluster_dic['SDSSJ1209+2640']           = {'ra':182.350000	,'dec':26.679722  ,'name':'SDSSJ1209+2640'}
    cluster_dic['SDSSJ0915:3826']           = {'ra':138.912500	,'dec':38.449722  ,'name':'SDSSJ0915:3826'}
    cluster_dic['SDSSJ0928+2031']           = {'ra':142.025000	,'dec':20.523889  ,'name':'SDSSJ0928+2031'}
    cluster_dic['SDSSJ0952+3434']           = {'ra':148.166667	,'dec':34.579722  ,'name':'SDSSJ0952+3434'}
    cluster_dic['SDSSJ0957+0509']           = {'ra':149.412500	,'dec':5.158889   ,'name':'SDSSJ0957+0509'}
    cluster_dic['SDSSJ1002+2031']           = {'ra':150.612500	,'dec':20.517500  ,'name':'SDSSJ1002+2031'}
    cluster_dic['SDSSJ1038+4849']           = {'ra':159.679167	,'dec':48.821944  ,'name':'SDSSJ1038+4849'}
    cluster_dic['SDSSJ1050+0017']           = {'ra':162.666667	,'dec':0.285278   ,'name':'SDSSJ1050+0017'}
    cluster_dic['SDSSJ1055+5547']           = {'ra':163.770833	,'dec':55.806389  ,'name':'SDSSJ1055+5547'}
    cluster_dic['SDSSJ1115+1645']           = {'ra':168.766667	,'dec':16.760556  ,'name':'SDSSJ1115+1645'}
    cluster_dic['SDSSJ1336-0331']           = {'ra':204.000000	,'dec':-3.524722  ,'name':'SDSSJ1336-0331'}
    cluster_dic['SDSSJ1138+2754']           = {'ra':174.537500	,'dec':27.908611  ,'name':'SDSSJ1138+2754'}
    cluster_dic['SDSSJ1152+0939']           = {'ra':178.195833	,'dec':9.504167   ,'name':'SDSSJ1152+0939'}
    cluster_dic['SDSSJ1152+3312']           = {'ra':178.000000	,'dec':33.228333  ,'name':'SDSSJ1152+3312'}
    cluster_dic['SDSSJ1156+1911']           = {'ra':179.025000	,'dec':19.186944  ,'name':'SDSSJ1156+1911'}
    cluster_dic['SDSSJ1207+5254']           = {'ra':181.900000	,'dec':52.916111  ,'name':'SDSSJ1207+5254'}
    cluster_dic['SDSSJ1632+3500']           = {'ra':248.041667	,'dec':35.008333  ,'name':'SDSSJ1632+3500'}
    cluster_dic['SDSSJ1522+2535']           = {'ra':230.720833	,'dec':25.594167  ,'name':'SDSSJ1522+2535'}
    cluster_dic['SDSSJ1329+2243']           = {'ra':202.391667	,'dec':22.716667  ,'name':'SDSSJ1329+2243'}
    cluster_dic['SDSSJ1343+4155']           = {'ra':205.891667	,'dec':41.910000  ,'name':'SDSSJ1343+4155'}
    cluster_dic['SDSSJ1420+3955']           = {'ra':215.162500	,'dec':39.918056  ,'name':'SDSSJ1420+3955'}
    cluster_dic['SDSSJ1439+1208']           = {'ra':219.791667	,'dec':12.140278  ,'name':'SDSSJ1439+1208'}
    cluster_dic['SDSSJ1456+5702']           = {'ra':224.004167	,'dec':57.039167  ,'name':'SDSSJ1456+5702'}
    cluster_dic['SDSSJ1527+0652']           = {'ra':231.937500	,'dec':6.875556   ,'name':'SDSSJ1527+0652'}
    cluster_dic['SDSSJ1531+3414']           = {'ra':232.795833	,'dec':34.240278  ,'name':'SDSSJ1531+3414'}
    cluster_dic['SDSSJ1621+0607']           = {'ra':245.383333	,'dec':6.122222   ,'name':'SDSSJ1621+0607'}
    cluster_dic['SDSSJ1723+3411']           = {'ra':260.900000	,'dec':34.199722  ,'name':'SDSSJ1723+3411'}
    cluster_dic['SDSSJ2111-0114']           = {'ra':317.829167	,'dec':-1.244444  ,'name':'SDSSJ2111-0114'}
    cluster_dic['SDSSJ2243-0935']           = {'ra':340.845833	,'dec':-9.589167  ,'name':'SDSSJ2243-0935'}
    cluster_dic['SDSSJ1604+2244']           = {'ra':241.041667	,'dec':22.738056  ,'name':'SDSSJ1604+2244'}


    Avalarr = np.zeros([len(cluster_dic.keys()),3])
    for kk, key in enumerate(cluster_dic.keys()):
        if verbose: print ' ---------------- Get extinction for '+cluster_dic[key]['name']+' ----------------'
        if verbose == 'Full':
            vb_getAv_area = True
        else:
            vb_getAv_area = False
        A, EBV, grid = kbs.getAv_area(cluster_dic[key]['ra'],cluster_dic[key]['dec'],radius,stepsize=stepsize,
                                      valreturn=valreturn,filter=filter,verbose=vb_getAv_area)
        Avalarr[kk,:] = A, EBV, grid
        if verbose: print ' - The results is A, EBV, grid = ',A,',', EBV,',', grid,','

    return Avalarr
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =