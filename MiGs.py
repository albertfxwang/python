"""
----------------------------
   NAME
----------------------------
 MiGs.py
----------------------------
   PURPOSE/DESCRIPTION
----------------------------
 MUSE inspection GUIs - GUIs for inspecting MUSE data products.
----------------------------
   COMMENTS
----------------------------

----------------------------
   EXAMPLES/USAGE
----------------------------


----------------------------
   BUGS
----------------------------

----------------------------
     FEATURES TO ADD
----------------------------

----------------------------
"""
#-------------------------------------------------------------------------------------------------------------
# IMPORTING MODULES
from Tkinter import *
import os
import sys
import glob
import datetime
import time
import numpy as np
import pdb
import subprocess
import pyfits
import scipy.ndimage
import commands
import MiGs
import collections
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
#-------------------------------------------------------------------------------------------------------------
def launch_MiG1D(directory='./',outputfile='DEFAULT',idsearchstr='spectrum_OBJID*.fits',
                 idlength=8,col_flux='FLUX',col_fluxerr='FLUXERR',col_wave='WAVE_AIR',
                 fluxunit='$10^{-20}$erg/s/cm$^2$/\\AA',lineuncertainty=False,
                 objlist=None,inspectorname='John Doe',clobber=False,infofile=None,col_infoid='UNIQUE_ID',
                 ds9xpa=False,openfitsauto=False,openfitsext='[0]',check4duplicates=False,skipempty=False,
                 outputcheck=False,latexplotlabel=False,autosaveplot=False,verbose=True,
                 skyspectrum='/Users/kschmidt/work/MUSE/skytable.fits'):
    """
    Launch the MUSE inspection GUI designed for looking at (any) 1D spectra following the MUSE-Wide format,
    i.e., binary fits table with columns wave_air, flux and fluxerr columns.
    """
    dir = directory
    if outputfile == 'DEFAULT':
        outfile = dir+'MiG1D_defaultoutput.txt'
    else:
        outfile = dir+outputfile
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # setup and launch GUI
    root = Tk()
    root.title("MUSE inspection GUI for 1D spectra (MiG1D)")
    root.geometry("910x430") # size of GUI window
    app = Application_1D(directory,outfile,master=root,idsearchstr=idsearchstr,idlength=8,col_flux='FLUX',
                         col_fluxerr=col_fluxerr,col_wave=col_wave,fluxunit=fluxunit,
                         objlist=objlist,verbose=verbose,ds9xpa=ds9xpa,openfitsauto=openfitsauto,openfitsext=openfitsext,
                         iname=inspectorname,latexplotlabel=latexplotlabel,clobber=clobber,infofile=infofile,
                         col_infoid=col_infoid,skipempty=skipempty,check4duplicates=check4duplicates,
                         outputcheck=outputcheck,autosaveplot=autosaveplot,skyspectrum=skyspectrum,
                         lineuncertainty=lineuncertainty)
    app.mainloop()
    root.destroy()
#-------------------------------------------------------------------------------------------------------------
def load_MiGoutput(MiGoutputfile,migversion='MiG1D',verbose=True):
    """
    Loading the output from MiG inspections

    """
    if migversion == 'MiG1D':
        if verbose: print ' - Loading \n   '+MiGoutputfile
        outputdata = np.genfromtxt(MiGoutputfile,comments='#',skip_header=2,names=True,dtype=None)

        comments = {}
        header   = []
        lines = open(MiGoutputfile,'r')
        for line in lines:
            line = line.replace('\r',' ')      # remove alternative end-of-line characters (return); only respect (\n)
            if line[0] == '#':                 # skipping header
                header.append(line[:-1])
            elif line == '\n':
                pass
            else:
                ID           = str(line.split()[0])
                comments[ID] = line.split('#C#')[-1][:-1]

    return outputdata, comments, header
#-------------------------------------------------------------------------------------------------------------
def getPID(searchstring,verbose=False):
    """
    Return PID of most recent process including the given search string
    Will ignore instances containing 'grep' and with a CPU time stamp of 0 seconds
    """
    cmd = "ps -eo pid,etime,command | grep "+searchstring
    fileobj = os.popen(cmd)            # return file object for ps command
    lines   = fileobj.readlines()      # read ouptu from ps command
    Nproces = len(lines)               # number of prcesses including searchstring

    PIDlist  = []
    time     = []
    for line in lines:
        if 'grep' not in line: # ignore the grep command PID
            ls       = line.split() # etime is [[dd-]hh:]mm:ss
            tsplit   = ls[1].split(':')
            if len(tsplit) == 2:    # if process has been running for minutes & seconds
                timeval = float(tsplit[0])*60. + float(tsplit[1])
            else:
                if '-' in tsplit[0]: # if process has been running for days
                    dayhour = tsplit[0].split('-')
                    timeval = float(dayhour[0])*24*60.*60 + \
                              float(dayhour[1])*60*60     + \
                              float(tsplit[1])*60. + float(tsplit[2])
                else:                # if process has been running for hours
                    timeval = float(tsplit[0])*60.*60     +\
                              float(tsplit[1])*60. + float(tsplit[2])

            if timeval > 0: # ignore 0.00 s instances
                if verbose: print 'Process:',line
                PIDlist.append(int(ls[0]))
                time.append(timeval)
            else:
                if verbose: print ' - Ignoring the following as it has a time stamp of 0s:'
                if verbose: print '   ',line
    if len(PIDlist) == 0:
        if verbose: print ' - No processes with given search string ('+searchstring+') in them. Returning None'
        return None

    if verbose: print 'PIDlist:',PIDlist
    if verbose: print 'time   :',time

    PID = np.array(PIDlist)[time == np.min(time)]
    if len(PID) > 1:
        print ' - Note multiple IDs with the same time stamp were found: ',PID
        print '   Returning the first PID'

    return PID[0]
#-------------------------------------------------------------------------------------------------------------
def check_idlist(idlist,dir,idsearchstr,verbose=True):
    """
    Checking if pngs exist for objects in idlist.
    Returning list of ids with existing files
    """
    if verbose: print ' - Checking ID list to make sure data for objects exists'
    goodids = np.array([])
    for objid in idlist:
        idstr = str(objid)
        pngs  = glob.glob(dir+'/'+idsearchstr.replace('OBJID',idstr))
        if len(pngs) > 0:
            goodids = np.append(goodids,objid)

    if (len(goodids) == 0):
        if verbose: print ' - WARNING None of the IDs have data in dir=\n   '+dir

    return goodids
#-------------------------------------------------------------------------------------------------------------
def get_objinfo(infofile,objid,idcol):
    """
    Return information on object given an input file
    """
    if infofile == None:
        returndat = None
    else:
        infodat = pyfits.open(infofile)[1].data
        objent  = np.where(infodat[idcol].astype(int) == int(objid))[0]

        if len(objent) == 0:
            returndat = None
        else:
            returndat = infodat[objent]

    return returndat
#-------------------------------------------------------------------------------------------------------------
def linelistdic(listversion='full'):
    """

    Dictionary containing lines with info on plotting etc.

    Line wavelength are taken from Morton 1991 when available

    """

    linelist = collections.OrderedDict()

    if listversion == 'full':
        #                         name                           wavelength[A]  horizontalalignment      lineref
        linelist['ovi1']   = ['OVI $\\lambda$1032'                , 1031.9261,         'right'      , 'Morton1991tab2']
        linelist['ovi2']   = ['OVI $\\lambda$1038'                , 1037.6167,         'left'       , 'Morton1991tab2']
        linelist['lya']    = ['Ly$\\alpha$ $\\lambda$1216'        , 1215.6737,         'right'      , 'Morton1991tab5']
        linelist['lyb']    = ['Ly$\\beta$ $\\lambda$1025'         , 1025.7219,         'right'      , 'Morton1991tab5']
        linelist['lyg']    = ['Ly$\gamma$ $\\lambda$973'          ,  972.5371,         'right'      , 'Morton1991tab5']
        linelist['nv1']    = ['NV $\\lambda$1239'                 , 1238.821 ,         'right'      , 'Morton1991tab5']
        linelist['nv2']    = ['NV $\\lambda$1243'                 , 1242.804 ,         'left'       , 'Morton1991tab5']
        linelist['cii']    = ['CII $\\lambda$1336'                , 1335.6627,         'right'      , 'Morton1991tab5']
        linelist['Siiv1']  = ['SiIV $\\lambda$1394'               , 1393.755 ,         'right'      , 'Morton1991tab5']
        linelist['Siiv2']  = ['SiIV $\\lambda$1403'               , 1402.770 ,         'left'       , 'Morton1991tab5']
        linelist['oiv1']   = ['OIV $\\lambda$1397'                , 1397.232 ,         'right'      , 'Morton1991tab5']
        linelist['oiv2']   = ['OIV $\\lambda$1400'                , 1399.780 ,         'left'       , 'Morton1991tab5']
        linelist['civ1']   = ['CIV $\\lambda$1548'                , 1548.195 ,         'right'      , 'Morton1991tab5']
        linelist['civ2']   = ['CIV $\\lambda$1551'                , 1550.770 ,         'left'       , 'Morton1991tab5']
        linelist['heii']   = ['HeII $\\lambda$1640'               , 1640.420 ,         'right'      , 'vandenberk+2001']
        linelist['oiiib1'] = ['OIII] $\\lambda$1661'              , 1660.809 ,         'right'      , 'Morton1991tab2']
        linelist['oiiib2'] = ['OIII] $\\lambda$1666'              , 1666.150 ,         'left'       , 'Morton1991tab2']
        linelist['ciii1']  = ['[CIII] $\\lambda$1907'             , 1907.    ,         'right'      , 'stark+2015']
        linelist['ciii2']  = ['CIII] $\\lambda$1909'              , 1909.    ,         'left'       , 'stark+2015']
        linelist['ciib']   = ['CII] $\\lambda$2326'               , 2326.113 ,         'right'      , 'Morton1991tab5']
        linelist['mgii1']  = ['MgII] $\\lambda$2796'              , 2795.528 ,         'right'      , 'Morton1991tab5']
        linelist['mgii2']  = ['MgII] $\\lambda$2803'              , 2802.705 ,         'left'       , 'Morton1991tab5']
        linelist['oii1']   = ['[OII] $\\lambda$3726'              , 3726.    ,         'right'      , 'Pradhan2006']
        linelist['oii2']   = ['[OII] $\\lambda$3729'              , 3729.    ,         'left'       , 'Pradhan2006']
        linelist['hd']     = ['H$\delta$ $\\lambda$4103'          , 4102.89  ,         'right'      , 'VandenBerk2001tab4']

        # Lines from http://www.sdss.org/dr7/algorithms/linestable.html,
        # http://adsabs.harvard.edu/abs/2008ApJS..174..282L, and
        # http://adsabs.harvard.edu/cgi-bin/nph-data_query?bibcode=2001AJ....122..549V&link_type=ABSTRACT
        #linelist['hd']     = ['$H\delta$'                         , 4101.74  ,         'right'      , '???']
        linelist['hg']     = ['H$\gamma$ $\\lambda$4340'          , 4340.47  ,         'right'      , '???']
        linelist['hb']     = ['H$\\beta$ $\\lambda$4861'          , 4861.33  ,         'right'      , '???']
        linelist['oiii1']  = ['[OIII] $\\lambda$4959'             , 4959.    ,         'right'      , '???']
        linelist['oiii2']  = ['[OIII] $\\lambda$5007'             , 5007.    ,         'right'      , '???']
        linelist['hei']    = ['HeI $\\lambda$5877'                , 5877.    ,         'right'      , '???']
        linelist['oi']     = ['OI $\\lambda$6302'                 , 6302.    ,         'right'      , '???']
        linelist['nii1']   = ['NII $\\lambda$6548'                , 6548.    ,         'right'      , '???']
        linelist['ha']     = ['H$\\alpha$ $\\lambda$6563'         , 6562.8   ,         'right'      , '???']
        linelist['nii2']   = ['NII $\\lambda$6583.5 '             , 6583.5   ,         'right'      , '???']
        linelist['sii1']   = ['SII $\\lambda$6718,'               , 6718.    ,         'right'      , '???']
        linelist['sii2']   = ['SII $\\lambda$6732'                , 6732.    ,         'right'      , '???']
        linelist['siii1']  = ['[SIII] $\\lambda$9071.1'           , 9071.1   ,         'right'      , '???']
        linelist['siii2']  = ['[SIII] $\\lambda$9533.2'           , 9533.2   ,         'right'      , '???']
    elif listversion == 'rest_uv_main':
        #                         name                           wavelength[A]  horizontalalignment      lineref
        # linelist['lyg']    = ['Ly$\gamma$ $\\lambda$973'          ,  972.5371,         'right'      , 'Morton1991tab5']
        linelist['ovi1']   = ['OVI $\\lambda$1032'                , 1031.9261,         'right'      , 'Morton1991tab2']
        linelist['ovi2']   = ['OVI $\\lambda$1038'                , 1037.6167,         'left'       , 'Morton1991tab2']
        linelist['lyb']    = ['Ly$\\beta$ $\\lambda$1025'         , 1025.7219,         'right'      , 'Morton1991tab5']
        linelist['lya']    = ['Ly$\\alpha$ $\\lambda$1216'        , 1215.6737,         'right'      , 'Morton1991tab5']
        linelist['nv1']    = ['NV $\\lambda$1239'                 , 1238.821 ,         'right'      , 'Morton1991tab5']
        linelist['nv2']    = ['NV $\\lambda$1243'                 , 1242.804 ,         'left'       , 'Morton1991tab5']
        # linelist['cii']    = ['CII $\\lambda$1336'                , 1335.6627,         'right'      , 'Morton1991tab5']
        # linelist['Siiv1']  = ['SiIV $\\lambda$1394'               , 1393.755 ,         'right'      , 'Morton1991tab5']
        # linelist['oiv1']   = ['OIV $\\lambda$1397'                , 1397.232 ,         'right'      , 'Morton1991tab5']
        # linelist['oiv2']   = ['OIV $\\lambda$1400'                , 1399.780 ,         'left'       , 'Morton1991tab5']
        # linelist['Siiv2']  = ['SiIV $\\lambda$1403'               , 1402.770 ,         'left'       , 'Morton1991tab5']
        linelist['civ1']   = ['CIV $\\lambda$1548'                , 1548.195 ,         'right'      , 'Morton1991tab5']
        linelist['civ2']   = ['CIV $\\lambda$1551'                , 1550.770 ,         'left'       , 'Morton1991tab5']
        # linelist['heii']   = ['HeII $\\lambda$1640'               , 1640.420 ,         'right'      , 'vandenberk+2001']
        # linelist['oiiib1'] = ['OIII] $\\lambda$1661'              , 1660.809 ,         'right'      , 'Morton1991tab2']
        # linelist['oiiib2'] = ['OIII] $\\lambda$1666'              , 1666.150 ,         'left'       , 'Morton1991tab2']
        linelist['ciii1']  = ['[CIII] $\\lambda$1907'             , 1907.    ,         'right'      , 'stark+2015']
        linelist['ciii2']  = ['CIII] $\\lambda$1909'              , 1909.    ,         'left'       , 'stark+2015']
        # linelist['ciib']   = ['CII] $\\lambda$2326'               , 2326.113 ,         'right'      , 'Morton1991tab5']
        # linelist['mgii1']  = ['MgII] $\\lambda$2796'              , 2795.528 ,         'right'      , 'Morton1991tab5']
        # linelist['mgii2']  = ['MgII] $\\lambda$2803'              , 2802.705 ,         'left'       , 'Morton1991tab5']
        # linelist['oii1']   = ['[OII] $\\lambda$3726'              , 3726.    ,         'right'      , 'Pradhan2006']
        # linelist['oii2']   = ['[OII] $\\lambda$3729'              , 3729.    ,         'left'       , 'Pradhan2006']
    else:
        sys.exit('invalid "listversion" provided ('+listversion+')')

    return linelist
#-------------------------------------------------------------------------------------------------------------
class Application_1D(Frame):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__(self,dir,outfile,master=None,idsearchstr='spectrum_OBJID*.fits',idlength=8,col_flux='FLUX',
                 col_fluxerr='FLUXERR',col_wave='WAVE_AIR',fluxunit='$10^{-20}$erg/s/cm$^2$/\\AA',
                 objlist=None,verbose=True,ds9xpa=False,openfitsauto=False,openfitsext='[0]',
                 iname='John Doe',latexplotlabel=False,clobber=False,infofile=None, col_infoid='UNIQUE_ID',
                 skipempty=False,check4duplicates=False,outputcheck=False,
                 autosaveplot=False,skyspectrum=False,lineuncertainty=False):
        """
        Intitialize the GUI

        -- INPUT --
        dir               Direcotory containing the data of the objects to inspect.
        outfile           Name of output file to create if it doesn't exists. Use clobber to overwrite.
        master            Provide 'master' display
        idsearchstr       Search string with 'OBJID' indicating the name of the ID in the filenames of the
                          1D spectra to look at in 'dir'. Used to generate 'objlist' if none provided.
        idlength          Length of ID used when extracting id's automatically from glob results on idserachstr
                          In this case, idsearchstr cannot contain a wildcard ('*') before the OBJID
        col_flux          Column name for fluxes in 1D spectrum
        col_fluxerr       Column name for flux errors in 1D spectrum
        col_wave          Column name for wavelengths in 1D spectrum
        fluxunit          String containing the flux unit (used to display on plot y-axis
        objlist           List of objects to inspect. If 'None' all objects in 'dir' will be
                          inspected.
        verbose           Toggle verbosity.
        ds9xpa            If xpa is availbale for comunicating commands to ds9
                          set this keyword to True and this will be used instead
                          of opening ds9 everytime the fits files are requested.

                          NB! XPA fixes the number of frames. If more than Nframes images are available they
                              will not be shown.

        openfitsauto      Automatically load the fits files into the DS9 window
                          when advancing to next (or previous) object.
        openfitsext       The fits extention to open when opening the fits files
        iname             Name of inspector to write in output file.
        latexplotlabel    Render plotting lables with latex; requires latex compiler.
        infofile          Fits catalog containing information on objects (linked by info_idcol)
        info_idcol        ID column in fits catalog provided in infofile
        clobber           Overwrites the output file if it already exists
        skipempty         Set to True to ignore unedited objects when writing to output file.
                          Hence, if skipempty = True objects with no comments, flags set or sliders changed
                          will be written to the output
        check4duplicates  Loop through output file whenever an object is save to check for
                          and remove duplicate entries
        outputcheck       Checking the written output to see if it contains the expected number
                          of objects etc.
        autosaveplot      Saving of the 1Dspec plot automatically when advancing to next object
        skyspectrum       Name of skyspectrum to plot (fits file with column flux - can be generated
                          with, e.g., http://www.eso.org/observing/etc/bin/simu/skycalc)
                          For wave < 10000A the MUSE sky will be plotted if it exists in "self.MUSEskydatdir"
        lineuncertainty   To plot shaded region around emission line indicators, to symbolize line position
                          uncertainty from redshift or line velocity offsets provide either Delta z or Delta v [km/s]
                          as uncertainty (lineuncertainty < 1 is treated as Delta z and lineuncertainty > 1 as Delta v)
        """
        self.now           = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.vb            = verbose
        self.idsearchstr   = idsearchstr
        self.col_flux      = col_flux
        self.col_fluxerr   = col_fluxerr
        self.col_wave      = col_wave
        self.fluxunit      = fluxunit
        self.dir           = dir
        self.infofile      = infofile
        self.col_infoid    = col_infoid
        self.quitting      = False
        self.duplicates    = check4duplicates
        self.outcheck      = outputcheck
        self.latex         = latexplotlabel
        self.autosaveplot  = autosaveplot
        self.skipempty     = skipempty
        self.skyspectrum   = skyspectrum
        self.MUSEskydatdir = '/Users/kschmidt/work/MUSE/skyspectra/'
        self.lineuncertainty = lineuncertainty

        self.ds9open       = False # set ds9 indicator (used for ds9xpa = False)
        self.ds9windowopen = False # set ds9 indicator (used for ds9xpa = True)
        self.xpa           = ds9xpa # check if user indacetd that xpa was available for ds9
        self.fitsauto      = openfitsauto # Open fits files automatically?
        self.openfitsext   = openfitsext
        if self.xpa:
            #sys.exit(' - XPA DS9 controls not enabled yet; still under construction (use ds9xpa=False)')
            self.ds9windowopen = False

        if skyspectrum:
            self.skydat = pyfits.open(skyspectrum)[1].data

        if os.path.exists(self.dir):
            self.twodfits = glob.glob(self.dir)
        else:
            sys.exit(' - The directory '+self.dir+' does not exist --> ABORTING')

        # -------- GET OBJIDS --------
        if objlist == None:
            self.file_1Dspec = [f for f in glob.glob(self.dir+'/'+idsearchstr.replace('OBJID','*'))]
            name_preID       = idsearchstr.split('OBJID')[0]
            self.objlist     = np.asarray([int(self.file_1Dspec[jj].split(name_preID)[-1][:idlength])
                                          for jj in xrange(len(self.file_1Dspec))])
        else:
            if type(objlist) == str:
                self.objlist = np.genfromtxt(objlist,dtype=None,comments='#')
            else:
                self.objlist = np.asarray(objlist)

            self.objlist = MiGs.check_idlist(self.objlist,self.dir,self.idsearchstr,verbose=self.vb)

        self.objlist    = np.unique(self.objlist).astype(int)

        if len(self.objlist) == 0:
            sys.exit(' No valid IDs found \n             Did you provide the right "idsearchstr"? ')

        self.currentobj = self.objlist[0]                    # set the first id to look at
        if verbose: print " - Found "+str(len(self.objlist))+' objects to inspect'

        # -------- OPEN/PREPARE OUTPUT FILE --------
        if os.path.isfile(outfile) & (clobber == True): # check if file is to be overwritten
            overwrite = raw_input(' - clobber==True Are you sure you want to overwrite '+outfile+'? (y/n): ')
            if (overwrite == 'y') or (overwrite == 'yes'):
                print "   Okay, I'll remove the file and start a new one"
                os.remove(outfile)
            elif (overwrite == 'n') or (overwrite == 'no'):
                print "   Okay, I'll append to the existing file, then"
            else:
                sys.exit('   "'+overwrite+'" is not a valid answer --> Aborting')

        if os.path.isfile(outfile):
            newfile   = False
            self.fout = open(outfile,'r')                # open existing file
            IDinspected = np.array([])                   # array to contain IDs in file
            for line in self.fout.readlines():           # loop through file to last line
                lsplit = line.split()
                if lsplit[0] != '#':
                    IDinspected = np.append(IDinspected,float(lsplit[0]))
            if len(IDinspected) == 0:
                sys.exit('Found no inspected objects in '+outfile)
            lastline = line
            self.fout.close()

            lastID = lastline.split()[0]                     # get the last ID in file
            if lastID != '#':
                objent = np.where(self.objlist == float(lastID))[0]
                if self.vb: print ' - The file '+outfile+' already exists (Resuming after last objects in output)'
                try:
                    self.currentobj = self.objlist[objent+1][0]  # change first id to look at
                except:
                    sys.exit(' - The last object in the outputfile is the last in "objlist" --> ABORTING ')
                Nremaining = len(self.objlist[objent+1:])
                Ninspected = len(np.unique(np.sort(IDinspected)))
                if self.vb:
                    print ' - Info from existing output: '
                    print '   '+str(Nremaining)+' of '+str(len(self.objlist))+' IDs still need to be expected'
                    print '   Found '+str(Ninspected)+' IDs already inspected in file'

            else:
                if self.vb: print ' - The file '+outfile+' already exists (append as last row does not contain ID)'
            self.fout     = open(outfile,'a')
        else:
            if self.vb: print ' - The file '+outfile+' was created (did not exist)'
            self.fout     = open(outfile,'w')
            self.fout.write('# Results from Visual Inspection of zfits initiated on '+self.now+' \n')
            self.fout.write('# Inspector: '+iname+' \n')
            newfile = True

        self.outfile = outfile

        # -------- ADD LABEL --------
        self.openpngs() # open pngs for first object
        position = [0,0,1]
        self.labelvar = StringVar()
        label = Label(master,textvariable=self.labelvar)
        label.grid(row=position[0],column=position[1],columnspan=position[2],sticky=N)
        self.labelvar.set(self.infostring())

        # -------- CREATE WIDGETS --------
        Frame.__init__(self, master)
        self.grid()
        self.create_widgets()

        # -------- SETUP DATAPLOT --------
        self.dataPlot_init(xsize=1200,ysize=600)
        self.dataPlot_loaddata()
        self.dataPlot_plot(refresh=False,newobj=True)
        self.DPxlow_full, self.DPxhigh_full, self.DPylow_full, self.DPyhigh_full = \
            self.dataPlot_getwindowinfo() # store first full window

        # -------- ADD IMAGE WINDOW --------
        # self.imgx,self.imgy = 990, 200
        # img = ImageTk.PhotoImage(Image.open(self.GUIimage).resize((self.imgx,self.imgy),Image.ANTIALIAS))
        # self.imageframe = Label(master, image=img)
        # self.imageframe.image = img
        # self.imageframe.grid(row = 150, column = 0, columnspan = 1, sticky=S)

        # -------- DRAW SEPERATORS --------
        self.drawsep(900,4,1 ,0,4,0,2,899,4)
        self.drawsep(900,4,29,0,4,0,2,899,4)
        self.drawsep(900,4,40,0,4,0,2,899,4)
        self.drawsep(900,4,60,0,4,0,2,899,4)
        self.drawsep(900,4,80,0,4,0,2,899,4)

        # -------- OPEN FITS FILES FOR FIRST OBJ --------
        if self.fitsauto: # loading fits files automatically
            if self.xpa:
                self.openfits_but_cmd_xpa()
            else:
                self.openfits_but_cmd()

        # -------- FINALIZE --------
        filehdr = '  '.join([key[3:] for key in self.keys])      # create header for output
        if newfile: self.fout.write('# ID '+filehdr+' byhandredshift byhandredshift_quality '
                                                       'multiple_redshift_solutions  \n')

        self.master.bind("<Key>", self.keyboard_cmd) # enable keyboard shortcuts
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def infostring(self):
        """
        Return string with information to display in GUI window
        """
        self.zobj = None # resetting redshift value
        objinfo   = MiGs.get_objinfo(self.infofile,self.currentobj,self.col_infoid)

        stringadd = ''
        if objinfo != None:
            cols      = objinfo.columns.names
            if 'FIELD_ID' in cols:
                fid = objinfo['FIELD_ID'][0]
                stringadd = stringadd+' field = '+str(fid)+'; '

            if 'IDENT_STRONGEST' in cols:
                line = objinfo['IDENT_STRONGEST'][0]
                stringadd = stringadd+' Strongest line = '+str(line)+'; '

            if 'REDSHIFT' in cols:
                self.zobj = objinfo['REDSHIFT'][0]
                zobjerr   = objinfo['REDSHIFT_ERR'][0]
                stringadd = stringadd+' z = '+str("%.6f" % self.zobj)+'+/- '+str("%.6f" % zobjerr)+'; '
        else:
            stringadd = stringadd+' ___NO INFOFILE PROVIDED___ '
        infostr = ">>> ID = "+str(self.currentobj)+'; '+stringadd+' <<<'

        return infostr
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def create_widgets(self):
        """
        Arrange the individual parts of the GUI
        postions are given as [row,column,span]
        """

        self.cbpos = [5,0,1]
        self.checkboxes(self.cbpos)
        self.commentfield([self.cbpos[0]+5,2,1])

        position = [35,0,3]
        textdisp = " Slider Values:     0) None;     1) tentative;     2) low-S/N;     3) high-S/N"
        label    = StringVar()
        txtlab   = Label(self,textvariable=label)
        label.set(textdisp)
        txtlab.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)

        self.openfits_but([65,3,1])

        self.prev_but([70,0,1])
        self.quit_but([70,1,1])
        self.skip_but([70,2,1])
        self.next_but([70,3,1])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def drawsep(self,width,height,row,col,colspan,xleftbottom,yleftbottom,xrighttop,yrighttop):
        """
        Draw a seperator
        """
        cv = Canvas(self, width=width, height=height)
        cv.grid(row = row, column = col, columnspan = colspan, sticky=N)
        cv.create_rectangle(xleftbottom, yleftbottom, xrighttop, yrighttop,fill='black')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def checkboxes(self,position,disable=False):
        """
        Checkboxes for keywords to assign to object
        """
        self.Ncol  = 4.

        self.sliders      = ['a','b','c','d','e','f','g','h']
        self.empty        = ['m','n','o','p',
                             #'c','d',
                             #'i','j'#,'k','l',
                             'q','r','s','t',
                             'u','v','w','x',
                             'y','z']
        self.calculations = []#['c','C','g','G','d','D','h','H','p','P','l','L']
        colors            = self.getcolors()

        # Note that letters in () enables sorting of boxes
        self.keys = {}
        self.keys['(a) Lyb_1025']   = 0
        self.keys['(b) OVI_1032' ]  = 0
        self.keys['(c) NV_1240']   = 0
        self.keys['(d) SiIVOIV_1399']  = 0
        #
        self.keys['(e) CIV_1549']  = 0
        self.keys['(f) HeII_1640']  = 0
        self.keys['(g) CIII_1909']  = 0
        self.keys['(h) MgII_2800']  = 0
        #
        self.keys['(i) close_counterpart']  = 0
        self.keys['(j) poor_match']  = 0
        #self.keys['(k) empty1']  = 0
        #self.keys['(l) empty2']  = 0

        #
        #self.keys['(m) '+self.gris1+'Mild_Contamination'] = 0
        #self.keys['(n) '+self.gris2+'Mild_Contamination'] = 0
        #self.keys['(o) '+self.gris2+'Contam_Defect']  = 0
        #self.keys['(p) Contamination_Level_Type']  = 0
        #
        #self.keys['(q) '+self.gris1+'Moderate_Contamination'] = 0
        #self.keys['(r) '+self.gris2+'Moderate_Contamination'] = 0
        #self.keys['(s) '+self.dirstr+'Defect']  = 0
        #self.keys['(t) empty7']  = 0
        #
        #self.keys['(u) '+self.gris1+'Severe_Contamination'] = 0
        #self.keys['(v) '+self.gris2+'Severe_Contamination'] = 0
        #self.keys['(w) '+self.dirstr+'Star']  = 0
        #self.keys['(x) I_have_no_idea']  = 0
        #
        #self.keys['(y) '+self.gris1+'Continuum'] = 0
        #self.keys['(z) '+self.gris2+'Continuum'] = 0

        if (sys.version_info[0] == 2) & (sys.version_info[1] == 7): # sort dictionary if running python 2.7
            import collections
            self.keys = collections.OrderedDict(sorted(self.keys.items()))
        else:
            print 'WARNING Python version not 2.7 so not sorting dictionary of keywords(1)'

        Nkey = 0
        self.cbdic     = {}
        self.sliderdic = {}
        for key in self.keys:
            rowval = position[0]+int(np.floor(Nkey/self.Ncol))
            if ('(a)' not in key) & ('(b)' not in key) & ('(c)' not in key) & ('(d)' not in key):
                rowval = rowval+1 # adjusting row value after first set of sliders
                if ('(e)' not in key) & ('(f)' not in key) & ('(g)' not in key) & ('(h)' not in key):
                    rowval = rowval+2 # adjusting row value after second set of sliders
            colval = position[1]+int(np.round((Nkey/self.Ncol-np.floor((Nkey/self.Ncol)))*self.Ncol))

            self.keys[key] = Variable()

            if key[1] in self.sliders:
                self.slider = Scale(self, from_=0, to=3,label=key,variable = self.keys[key],
                                    orient=HORIZONTAL,background=colors[key[1]],length=200)
                self.slider.grid(row=rowval,column=colval,columnspan=position[2],rowspan=2,sticky=W)
                self.slider.set(0)

                if disable:
                    self.slider.configure(state='disabled')
                else:
                    self.sliderdic[key] = self.slider
            elif key[1] in self.empty:
                self.cb = Checkbutton(self, text='emptyXX')
                #self.cb.grid(row=position[0]+3,column=3,columnspan=1,sticky=W)
                self.cb.deselect()
                self.keys[key].set('-1')
                if key[1] in self.calculations:
                    self.keys[key].set(key)
            else:
                self.cb = Checkbutton(self, text=key, variable=self.keys[key],background=colors[key[1]])
                self.cb.grid(row=rowval,column=colval,columnspan=position[2],sticky=W)
                self.cb.deselect()

                if disable:
                    self.cb.configure(state='disabled')
                else:
                    self.cbdic[key] = self.cb

            Nkey = Nkey + 1
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def getcolors(self,):
        """
        Dictionary with colors for keys
        """
        collist = ['orange','red','cyan','magenta','green','white','gray']
        colors  = {}
        colors['a'] = collist[6]
        colors['b'] = collist[6]
        colors['c'] = collist[6]
        colors['d'] = collist[6]
        colors['e'] = collist[6]
        colors['f'] = collist[6]
        colors['g'] = collist[6]
        colors['h'] = collist[6]
        colors['i'] = collist[5]
        colors['j'] = collist[5]
        colors['k'] = collist[0]
        colors['l'] = collist[0]
        colors['m'] = collist[0]
        colors['n'] = collist[0]
        colors['o'] = collist[0]
        colors['p'] = collist[0]
        colors['q'] = collist[0]
        colors['r'] = collist[0]
        colors['s'] = collist[0]
        colors['t'] = collist[0]
        colors['u'] = collist[0]
        colors['v'] = collist[0]
        colors['w'] = collist[0]
        colors['x'] = collist[0]
        colors['y'] = collist[0]
        colors['z'] = collist[0]

        return colors

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_init(self,xsize=500,ysize=150,rowval=45):
        """
        Inititalize the data plot window (also used to kill it when exiting)
        """
        #----------------- Plot setup -----------------
        plt.ioff()                                                 # turn off interactive plotting
        plt.ion()                                                  # enable interactive plotting
        self.DPFsize  = 16
        self.DPlwidth = 2
        self.DPxscale = 1e4
        self.DPcolor  = ['blue','red','purple','magenta','salmon','orange']
        self.DPxrange = [0.4,1.7]
        if self.latex:
            plt.rc('text', usetex=True)                            # enabling LaTex rendering of text
        else:
            plt.rc('text', usetex=False)                           # disabling LaTex rendering of text
        plt.rc('font' , family='serif',size=self.DPFsize)          # setting text font
        plt.rc('xtick', labelsize=self.DPFsize)
        plt.rc('ytick', labelsize=self.DPFsize)

        self.dataPlot_fig    = plt.figure()
        self.dataPlot_fig.canvas.set_window_title('1D spectrum of object '+str(self.currentobj))
        self.dataPlot_fig.subplots_adjust(wspace=0.2, hspace=0.2,left=0.1, right=0.98, bottom=0.15, top=0.95)
        self.dataPlot_ax     = self.dataPlot_fig.add_subplot(111)
        self.dataPlotManager = plt.get_current_fig_manager()  # get plotting canvas
        self.dataPlotManager.resize(xsize,ysize)

        # ==== SLIDERS =====
        zinit = 0.0
        if self.zobj: zinit = self.zobj

        self.varsliderz = DoubleVar()
        self.sliderz  = Scale(self, from_=0.00, to=6.0,label='Redshift (n)+ (N)-',variable = self.varsliderz,
                              orient=HORIZONTAL,background='gray',length=200,resolution=0.001)
        self.sliderz.grid(row=rowval,column=0,columnspan=1,rowspan=1,sticky=W)
        self.varsliderz.set(zinit) # set intial value of slider

        self.varslidersmooth  = DoubleVar()
        self.slidersmooth= Scale(self, from_=0, to=15,label='Gauss smooth (m)+ (M)-',
                                 variable = self.varslidersmooth,
                                 orient=HORIZONTAL,background='gray',length=200,resolution=0.1)
        self.slidersmooth.grid(row=rowval,column=1,columnspan=1,rowspan=1,sticky=W)
        self.varslidersmooth.set(0) # set intial value of slider

        self.varsliderzqual = DoubleVar()
        self.sliderzqual    = Scale(self, from_=0, to=4.0,label='(q) By-hand redshift quality',
                                    variable = self.varsliderzqual,orient=HORIZONTAL,background='gray',
                                    length=200,resolution=1.0)
        self.sliderzqual.grid(row=rowval,column=2,columnspan=1,rowspan=1,sticky=W)
        self.varsliderzqual.set(0) # set intial value of slider


        # ==== COMMENT FIELD ====
        self.byhandzlabel = Label(self,text='(u) By-hand redshift:  ',background='white')
        self.byhandzlabel.grid(row=rowval,column=3,columnspan=1,sticky=NW)
        self.byhandz = Entry(self)
        self.byhandz.grid(row=rowval,column=3,columnspan=1,sticky=SW)

        # ==== CHECK BOX ====
        self.skyboxvar = Variable()
        self.skybox = Checkbutton(self, text='(o) Show Sky', variable=self.skyboxvar,background='white')
        self.skybox.grid(row=rowval+1,column=0,columnspan=1,sticky=W)
        self.skybox.deselect()

        self.err1Dboxvar = Variable()
        self.err1Dbox    = Checkbutton(self, text='(p) Show 1D errors', variable=self.err1Dboxvar,
                                       background='white')
        self.err1Dbox.grid(row=rowval+1,column=1,columnspan=1,sticky=W)
        self.err1Dbox.deselect()
        #if (self.GiGf == None): self.err1Dbox.configure(state='disabled')

        self.mzsboxvar = Variable()
        self.mzsbox    = Checkbutton(self, text='(t) Multiple Redshift Solutions', variable=self.mzsboxvar,
                                     background='white')
        self.mzsbox.grid(row=rowval+1,column=2,columnspan=1,sticky=W)
        self.mzsbox.deselect()

        # ==== BUTTONS ====
        self.dataPlot_fullzoombutton([rowval+2,0,1])
        self.dataPlot_redrawbutton([rowval+2,1,1])
        self.dataPlot_savebutton([rowval+2,3,1])

        self.DPxlow, self.DPxhigh, self.DPylow, self.DPyhigh = self.dataPlot_getwindowinfo() # store window
        # self.dataPlotManager.destroy()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_redrawbutton(self,position):
        """
        Button to redraw plot
        """
        self.dpbut_redraw = Button(self)
        self.dpbut_redraw["text"] = "(r) Redraw"
        self.dpbut_redraw["command"] = self.dataPlot_redrawbutton_cmd
        self.dpbut_redraw.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_redrawbutton_cmd(self):
        """
        Command for redrawing the plot
        """
        self.DPxlow, self.DPxhigh, self.DPylow, self.DPyhigh = self.dataPlot_getwindowinfo() # store window
        self.dataPlot_plot(refresh=True,verbose=True)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_fullzoombutton(self,position):
        """
        Button to go to full zoom in the plot
        """
        self.dpbut_fullzoom = Button(self)
        self.dpbut_fullzoom["text"] = "(z) full zoom"
        self.dpbut_fullzoom["command"] = self.dataPlot_fullzoombutton_cmd
        self.dpbut_fullzoom.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_fullzoombutton_cmd(self):
        """
        Command for going back to full zoom in the plot
        """
        self.DPxlow, self.DPxhigh, self.DPylow, self.DPyhigh = self.dataPlot_getwindowinfo() # store window
        self.dataPlot_plot(refresh=True,fullzoom=True,verbose=True)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_savebutton(self,position):
        """
        Button to save the plot using default naming
        """
        self.dpbut_save = Button(self)
        self.dpbut_save["text"] = "(s) Quick save plot"
        if self.autosaveplot: self.dpbut_save["text"] = "(s) Autosave Enabled"
        self.dpbut_save["command"] = self.dataPlot_savebutton_cmd
        self.dpbut_save.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_savebutton_cmd(self):
        """
        Command for saving the created plot
        """
        plotname = self.dir+str("%.5d" % self.currentobj)+'_MiG1D_specplot.pdf'
        self.dataPlot_fig.savefig(plotname)
        print ' - Saved plot window to \n   '+plotname
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_loaddata(self,verbose=False):
        """
        Loading the data for dataPlot so it's only necessary to load data once

        """
        self.DPidstr         = str("%.5d" % self.currentobj)
        globresult           = glob.glob(self.dir+'/'+self.idsearchstr.replace('OBJID',str(self.currentobj)))
        self.fits1Dfound     = [ff for ff in globresult if ('_2D' not in ff)] # ignore *_2D* fits
        self.fitsallfound    = glob.glob(self.dir+'/*'+str(self.currentobj)+'*.fits')
        self.DP_wave_all     = []
        self.DP_flux_all     = []
        self.DP_fluxerr_all  = []
        self.DP_contam_all   = []
        for f1D in self.fits1Dfound:
            dat1D   = pyfits.open(f1D)[1].data
            self.DP_wave_all.append(dat1D[self.col_wave])
            self.DP_flux_all.append(dat1D[self.col_flux])
            self.DP_fluxerr_all.append(dat1D[self.col_fluxerr])
            if 'CONTAM' in dat1D.columns.names:
                self.DP_contam_all.append(dat1D['CONTAM'])
            else:
                self.DP_contam_all.append(dat1D[self.col_flux]*0.0-99.)

        objinfo = MiGs.get_objinfo(self.infofile,self.currentobj,self.col_infoid)
        if objinfo != None:
            cols      = objinfo.columns.names
            if 'REDSHIFT' in cols:
                self.varsliderz.set(objinfo['REDSHIFT'][0])

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    def dataPlot_plot(self,verbose=False,refresh=False,newobj=False,fullzoom=False):
        """
        Plotting the 1D spectra loaded in dataPlot_loaddata

        """
        self.dataPlot_fig.canvas.set_window_title('1D spectrum of object '+str(self.currentobj))
        xlow, xhigh, ylow, yhigh = self.dataPlot_getwindowinfo()
        if fullzoom:
            xlow, xhigh, ylow, yhigh =  self.DPxlow_full, self.DPxhigh_full, self.DPylow_full, self.DPyhigh_full
        #----------------- Define emission line list -----------------
        llist    = MiGs.linelistdic(listversion='full')
        linelist = np.asarray([llist[key][1] for key in llist.keys()])
        linename = [llist[key][0] for key in llist.keys()]
        #----------------- Refreshing plot window-----------------
        if refresh:
            self.dataPlot_fig.clf() # clearing figure
            self.dataPlot_ax     = self.dataPlot_fig.add_subplot(111)
        #----------------- Grab info from sliders -----------------
        smoothlevel  = float(self.varslidersmooth.get())
        if verbose: print ' - Grabbed the Gauss smooth level ',smoothlevel,' from the slider'
        redshift     = float(self.varsliderz.get())
        if verbose: print ' - Grabbed the redshift '+str("%.3f" % redshift)+' from the slider'

        try:
            zbyhand      = float(self.byhandz.get())
            if type(zbyhand) == float:
                redshift = zbyhand
                if verbose: print '   But the redshift',zbyhand,'was found in "by-hand" field so using that instead '
            self.varsliderz.set(zbyhand)
        except:
            pass

        #----------------- Flambda spec -----------------
        xrangeflam = self.DPxrange
        contamplotted = False
        ymax   = []
        ymin   = []

        for ii in range(len(self.fits1Dfound)):
            color     = self.DPcolor[ii]
            wave1D    = self.DP_wave_all[ii]/self.DPxscale # wavelengths converted from A to micron
            flux1D    = self.DP_flux_all[ii]
            flux1Derr = self.DP_fluxerr_all[ii]
            contam    = self.DP_contam_all[ii]

            if (len(flux1D) >= 1):
                ymin.append(np.min(flux1D))
                ymax.append(np.max(flux1D))
                labstr = self.fits1Dfound[ii].split('/')[-1]
                if self.latex:
                    labstr = labstr.replace('_','\_')
                # - - - - - - - - - - Spectrum itself - - - - - - - - - -
                self.dataPlot_ax.plot(wave1D, flux1D, color=color,linestyle='-',
                                      linewidth=self.DPlwidth*1.5, alpha=0.40)

                # - - - - - - - - - - Smoothed spectrum - - - - - - - - - -
                filtersigma   = smoothlevel
                flux1D_smooth = scipy.ndimage.filters.gaussian_filter1d(flux1D, filtersigma,cval=0.0)
                self.dataPlot_ax.plot(wave1D, flux1D_smooth, color=color,linestyle='-',
                                      linewidth=self.DPlwidth*1.5, alpha=1.0,label=labstr)

                # - - - - - Shaded error region around curve if requested - - - - -
                if (self.err1Dboxvar.get() != '0'):
                    xwinmin, xwinmax, ywinmin, ywinmax = self.dataPlot_getwindowinfo()
                    serr     = flux1Derr
                    serr[flux1Derr > 1e3] = 1e3 # fix errors to prevent "OverflowError: Allocated too many blocks"
                    filllow  = np.clip(flux1D,ywinmin,ywinmax)-serr
                    fillhigh = np.clip(flux1D,ywinmin,ywinmax)+serr
                    plt.fill_between(wave1D,filllow,fillhigh,alpha=0.20,color=color)

                # - - - - - - - - - - Contam curve is present - - - - - - - - - -
                if (contam != -99).any():
                    self.dataPlot_ax.plot(wave1D, contam, color=color,linestyle='--',
                                          linewidth=self.DPlwidth, alpha=0.40)
                    contam_smooth = scipy.ndimage.filters.gaussian_filter1d(contam, filtersigma,cval=0.0)
                    self.dataPlot_ax.plot(wave1D, contam_smooth, color=color,linestyle='--',
                                          linewidth=self.DPlwidth, alpha=1.0)
                    contamplotted = True

                # - - - - - - - - - - Sky spectrum - - - - - - - - - -
                if (self.skyboxvar.get() != '0'):
                    objinfo  = MiGs.get_objinfo(self.infofile,self.currentobj,self.col_infoid)
                    if (np.max(wave1D) < 1.0) & (len(objinfo) == 1):
                        fieldno          = objinfo['FIELD_ID']
                        skyMUSEfilename  = glob.glob(self.MUSEskydatdir+'SKY*cdfs*-'+str("%.2d" % fieldno)+'*av.fits')
                        skyMUSE          = pyfits.open(skyMUSEfilename[0])[1].data
                        skywave = skyMUSE['lambda']/self.DPxscale
                        skyent  = np.where((skywave > np.min(wave1D)) & (skywave < np.max(wave1D)))[0]
                        skywave = skywave[skyent]
                        skylow  = np.zeros(len(skywave))
                        skyflux = skyMUSE['data'][skyent]
                        skyhigh = skyflux/1.0
                    elif self.skyspectrum:
                        skywave = self.skydat['lam']
                        skyent  = np.where((skywave > np.min(wave1D)) & (skywave < np.max(wave1D)))[0]
                        skywave = skywave[skyent]
                        skylow  = np.zeros(len(skywave))
                        skyflux = self.skydat['flux'][skyent]
                        skymax  = np.sort(flux1D)[np.round(len(flux1D)*0.95)]
                        skyhigh = skyflux/np.max(skyflux)*skymax
                    else:
                        skywave = wave1D
                        skylow  = np.zeros(len(skywave))-100.
                        skyhigh = np.zeros(len(skywave))+100.

                    plt.fill_between(skywave,skylow,skyhigh,alpha=0.3,color='black')
                    skyhigh_smooth = scipy.ndimage.filters.gaussian_filter1d(skyhigh, filtersigma,cval=0.0)
                    plt.fill_between(skywave,skylow,skyhigh_smooth,alpha=0.8,color='black')

        # set ranges based on spectra
        if (len(ymin) != 0) & (len(ymax) != 0):
            yrangeflam = [0.95*min(ymin), 1.05*max(ymax)]
            # if yrangeflam[0] < -0.01: yrangeflam[0] = -0.01
            # if yrangeflam[1] >  10.0: yrangeflam[1] =  10.0
        else:
            yrangeflam = 0.0, 1000.0

        if not newobj: # only check window if not plotting new object
            if (ylow != yrangeflam[0]) or (yhigh != yrangeflam[1]):
                yrangeflam = [ylow,yhigh]
        Dyrange    = yrangeflam[1]-yrangeflam[0]
        self.dataPlot_ax.set_ylim(yrangeflam)

        if not newobj: # only check window if not plotting new object
            if (xlow != xrangeflam[0]) or (xhigh != xrangeflam[1]):
                xrangeflam = [xlow,xhigh]
        self.dataPlot_ax.set_xlim(xrangeflam)

        if self.latex:
            xlab = '$\lambda / [\mu\mathrm{m}]$'
            ylab = '$f_\\lambda$ / ['+self.fluxunit+']'
        else:
            xlab = 'lambda / [micron]'
            ylab = 'f_lambda / ['+self.fluxunit+']'

        self.dataPlot_ax.set_xlabel(xlab)
        self.dataPlot_ax.set_ylabel(ylab)

        self.dataPlotManager.canvas.draw()

        # === plot emission lines for scale ===
        for ii in range(len(linelist)):
            lineposition = linelist[ii]/self.DPxscale*(redshift+1.0)
            self.dataPlot_ax.plot(np.zeros(2)+lineposition,yrangeflam,color='#006600',alpha=0.7,
                                  linestyle='-',linewidth=self.DPlwidth*2)

            if self.lineuncertainty:
                if (self.lineuncertainty <= 1.0) & (self.lineuncertainty > 0.0):   # treat input as Delta z uncertainty
                    zoffset  = self.lineuncertainty
                elif self.lineuncertainty > 1.0: # treat input as Delta v uncertainty
                    zoffset = self.lineuncertainty*(redshift+1.0) / 299792.458
                else:
                    if self.vb: print ' WARNING: Invalid value of "lineuncertainty" using dz=0.1'
                    zoffset  = 0.1

                linexmin = ( (redshift-zoffset) +1) * linelist[ii]/self.DPxscale
                linexmax = ( (redshift+zoffset) +1) * linelist[ii]/self.DPxscale
                lineymin = yrangeflam[0]
                lineymax = yrangeflam[1]

                plt.fill_between(np.asarray([linexmin,linexmax]),np.zeros(2)+lineymin,np.zeros(2)+lineymax,
                                 alpha=0.2,color='#006600')

                voffset = zoffset * 299792.458 / (redshift+1.0)

            textpos = linelist[ii]/self.DPxscale*(redshift+1.0)
            if (textpos > xrangeflam[0]) & (textpos < xrangeflam[1]):
                self.dataPlot_ax.text(textpos,yrangeflam[0]+Dyrange*0.05,
                                      linename[ii],color='#006600',size=self.DPFsize-3.,rotation='vertical',
                                      horizontalalignment='right',verticalalignment='bottom')

        # === position legend ===
        box = self.dataPlot_ax.get_position()
        self.dataPlot_ax.set_position([box.x0, box.y0, box.width, box.height * 0.83])



        if (self.skyboxvar.get() != '0'):
            self.dataPlot_ax.plot(0,0,'black',alpha=0.8,label='Sky spectrum',linewidth=self.DPlwidth*2)
        if contamplotted:
            self.dataPlot_ax.plot(0,0,'black',alpha=0.8,label='Contamination',linewidth=self.DPlwidth,ls='--')
        self.dataPlot_ax.plot(0,0,'green',label='Lines at $z$ = '+str("%.3f" % redshift),linewidth=self.DPlwidth*2)
        if self.lineuncertainty:
            linelab = 'Line uncertainty $z$ +/- '+str("%.4f" % zoffset)+' (+/- '+str("%.f" % voffset)+' km/s)'
            self.dataPlot_ax.plot(0,0,'#006600',alpha=0.4,label=linelab,linewidth=self.DPlwidth*5)

        leg = self.dataPlot_ax.legend(fancybox=True, loc='upper center',numpoints=1,prop={'size':self.DPFsize-3.},
                                      ncol=2,bbox_to_anchor=(0.5, 1.27))
        #leg.get_frame().set_alpha(0.7)

        self.dataPlotManager.canvas.draw()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def dataPlot_getwindowinfo(self):
        """
        get information about window after zoom etc.
        """
        xmin  = self.dataPlot_ax.get_xbound()[0]
        xmax  = self.dataPlot_ax.get_xbound()[1]
        ymin  = self.dataPlot_ax.get_ybound()[0]
        ymax  = self.dataPlot_ax.get_ybound()[1]

        return xmin, xmax, ymin, ymax
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def commentfield(self,position):
        """
        Field to provide comments
        """
        self.label = Label(self,text='(x) Comments ("tab" to move focus):  ')
        self.label.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
        self.comments = Entry(self)
        self.comments.grid(row=position[0],column=position[1]+position[2],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def openpngs(self,objid=None):
        """
        Function to open pngs (and pdfs) of object
        """
        if objid == None:
            id = self.currentobj
        else:
            id = objid
        idstr     = str(id)
        self.pngs = glob.glob('/Users/kschmidt/work/MUSE/ELpostagestamps/'+idstr+'*.png')+\
                    glob.glob(self.dir+'/*'+idstr+'*.png')+glob.glob(self.dir+'/*'+idstr+'*.pdf')

        if len(self.pngs) == 0:
            if self.vb: print ' - Did not find any png or pdf files to open. Looked for '+\
                              self.dir+'/*'+idstr+'*.pdf/pdf'
        else:
            self.file = self.pngs[0].split('/')[-1]

            pngorderedlist = self.pngs

            self.plat = sys.platform
            if self.plat == 'darwin':
                import platform
                macversion = platform.mac_ver()[0]
                if float(macversion.split('.')[1]) > 6: # check if "open -F" is available (mac OS X 10.7.0 and above)
                    opencmd = 'open -n -F '+' '.join(pngorderedlist)
                else:
                    opencmd = 'open -n '+' '.join(pngorderedlist)
            elif self.plat == 'linux2' or 'Linux':
                opencmd = 'gthumb '+' '.join(pngorderedlist)+' &'
            else:
                sys.exit(' Did not recognize the platform you are running MiG1D on '
                         '(set up for "darwin" (MacOSX), linux and linux2) --> ABORTING')

            self.pPNG = subprocess.Popen(opencmd,shell=True,executable=os.environ["SHELL"])
            time.sleep(1.1)# sleep to make sure png appear in PIDlist
            if self.plat == 'darwin':
                self.pngPID = MiGs.getPID('Preview.app',verbose=False) # get PID of png process
            elif self.plat == 'linux2' or 'Linux':
                self.pngPID = MiGs.getPID('gthumb',verbose=False)      # get PID of png process
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def openfits_but(self,position):
        """
        Button to open fits files
        """
        self.regiontemp = 'temp_ds9_forinspection.reg'
        self.fitsb = Button(self)
        self.fitsb["text"] = "(0) Open fits files"
        if self.xpa:
            self.fitsb["command"] = self.openfits_but_cmd_xpa
        else:
            self.fitsb["command"] = self.openfits_but_cmd

        self.fitsb.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def openfits_but_cmd_xpa(self):
        """
        Command for openfits button
        """
        lockstr = self.lockds9string()
        ds9cmd  = ' '

        if not self.ds9windowopen:
            ds9cmd = ds9cmd+'ds9 -geometry 1200x600 -scale zscale '+\
                     lockstr+' -tile grid layout 4 2 '
            self.pds9   = subprocess.Popen(ds9cmd,shell=True,executable=os.environ["SHELL"])
            time.sleep(1.1)# sleep to make sure ds9 appear in PIDlist
            self.ds9PID = MiGs.getPID('ds9',verbose=False) # get PID of DS9 process
            self.ds9windowopen = True
            time.sleep(1.0)
            for ii in np.arange(1,17):
                out = commands.getoutput('xpaset -p ds9 frame new')
            out = commands.getoutput('xpaset -p ds9 tile')

        Fstart = 1
        for ff, filename in enumerate(self.fitsallfound):
            imgname    = filename.split('/')[-1].replace('.fits','')
            out        = commands.getoutput('xpaset -p ds9 frame '+str(Fstart))
            regionfile = self.regiontemp.replace('.reg',imgname+'file'+str(ff)+'.reg')
            self.ds9textregion(imgname,filename=regionfile)
            out = commands.getoutput('xpaset -p ds9 file '+filename+self.openfitsext)
            out = commands.getoutput('xpaset -p ds9 regions '+regionfile)
            Fstart += 1

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def openfits_but_cmd(self):
        """
        Command for openfits button
        """
        self.ds9open = True
        lockstr = self.lockds9string()

        ds9cmd  = 'ds9 -geometry 1200x600 -scale zscale '+lockstr+' -tile grid layout 4 2'

        for ff, filename in enumerate(self.fitsallfound):
            imgname    = filename.split('/')[-1].replace('.fits','')
            regionfile = self.regiontemp.replace('.reg',imgname+'file'+str(ff)+'.reg')
            self.ds9textregion(imgname,filename=regionfile)
            ds9cmd = ds9cmd+' "'+filename+self.openfitsext+'" -region '+regionfile+' '

        self.pds9   = subprocess.Popen(ds9cmd,shell=True,executable=os.environ["SHELL"])
        time.sleep(1.1)# sleep to make sure ds9 appear in PIDlist
        self.ds9PID = MiGs.getPID('ds9',verbose=False) # get PID of DS9 process
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def lockds9string(self):
        """
        """
        # if int(self.ds9version[1].split('.')[0]) >= 7: # only lock if ds9 version is 7 or later
        #     lockstr = ' -lock frame physical '
        # else:
        #     print ' - WARNING DS9 version older than 7.*; Not locking frames.'
        #     lockstr = ' '

        lockstr = ' -lock frame physical '

        return lockstr
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def ds9textregion(self,text,filename='temp.reg'):
        """
        Create ds9 region file with text string
        Note that it's overwriting any existing file!
        """
        regstr = 'physical\n# text(130,10) textangle=0 textrotate=0 font="helvetica 12 normal roman" text={'+text+'}'
        fds9region = open(filename,'w')
        fds9region.write(regstr)
        fds9region.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def reset(self,skip=False):
        """
        Writing results to output, resetting checkboxes, and closing DS9 and PNG windows

        if skip=True nothing will be written to output file.
        """
        if (self.autosaveplot) & (skip==False): self.dataPlot_savebutton_cmd() # saving plot before resetting

        try: # checking that the input can be converted to a float
            zbyhand      = str(float(self.byhandz.get()))+' '
        except:
            zbyhand      = '-99 '
            if (str(self.byhandz.get()) != ''):
                print ' - WARNING: by-hand redshift field ('+str(self.byhandz.get())+\
                      ') could not be converted to float. Using -99.'
        zbyhand = zbyhand+str(int(self.varsliderzqual.get()))

        resultstr  = ' '+str("%.5d" % self.currentobj)+' '
        defaultstr = resultstr
        for key in self.keys:
            keyval    = str(self.keys[key].get())
            if keyval == '-1':
                defaultstr = defaultstr+' '+str(keyval)
            elif len(keyval) > 10: # for text keys
                defaultstr = defaultstr+' '+keyval
            else:
                defaultstr = defaultstr+' '+str(0)
            resultstr = resultstr+' '+str(keyval)

        # by-hand redshift info
        defaultstr = defaultstr+' -99 0'
        resultstr  = resultstr+' '+zbyhand

        # Multiple redshift solutions?
        defaultstr = defaultstr+' 0'
        resultstr  = resultstr+' '+self.mzsboxvar.get()

        # adding info from comment and wave fields
        defaultstr = defaultstr +'  #C#  \n'
        resultstr  = resultstr  +'  #C# '+self.comments.get()+' \n'

        skipin = skip # storing original skip value
        if (resultstr == defaultstr) & (self.skipempty == True): skip = True
        if not skip:
            if self.duplicates:
                Ndup = self.removeoutputduplicate(self.currentobj)

            self.fout.write(str(resultstr))
        if resultstr == defaultstr: skip = skipin # restoring original skip value

        # --- close and re-open output file so inspection is saved ---
        self.fout.close()
        self.fout = open(self.outfile,'a')

        # --- resetting widgets and closing windows ---
        self.comments.delete(0,END) # reset comment field
        self.byhandz.delete(0,END)

        self.varsliderz.set(0.0) # set intial value of slider
        self.varslidersmooth.set(0)    # set intial value of slider
        self.varsliderzqual.set(0)     # set intial value of slider

        self.checkboxes(self.cbpos) # reset check boxes
        self.skybox.deselect()
        self.err1Dbox.deselect()
        self.mzsbox.deselect()

        self.closewindows()

        self.ds9open = False # resetting ds9 indicator
        self.focus_set() # set focus to main window
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def removeoutputduplicate(self,id):
        """
        Subtract continuum from science fram
        """
        self.fout.close()
        idstr       = str("%.5d" % id)
        stringstart = ' '+idstr+' '
        file        = open(self.outfile,'r')
        lines       = file.readlines()
        file.close()
        file = open(self.outfile,"w")

        Ndup        = 0
        for line in lines:
            if line[0:10] != stringstart:
                file.write(line)
            else:
                if self.vb: print ' - Found dublicate entry for ID '+idstr+' -> deleting it!'
                Ndup = Ndup+1

        file.close()
        self.fout = open(self.outfile,'a')
        return Ndup
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def closewindows(self):
        """
        Close PNG and DS9 windows
        """
        killsignal = 1     # see bash> man kill
        try:
            os.kill(self.pngPID,killsignal)                  # close PNG window for currentobj
        except:
            print '   WARNING error occurred while trying to close PNG window(s)'

        if np.logical_or(((self.ds9open == True) & (self.xpa == False)),
                         ((self.xpa == True) & (self.quitting == True) & (self.ds9windowopen == True))):
            try:
                os.kill(self.ds9PID,killsignal)                  # close DS9 window for currentobj
            except:
                if self.vb: print ' - WARNING: Could not kill DS9 process id ',self.ds9PID
            rmout = commands.getoutput('rm '+self.regiontemp.replace('.reg','*.reg')) # removing ds9 region file

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def skip_but(self,position):
        self.skip = Button(self)
        self.skip["text"] = "Skip object"
        self.skip["command"] = self.skip_but_cmd
        self.skip.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def skip_but_cmd(self):
        """
        Command for skip button
        """
        self.reset(skip=True)

        if self.currentobj == self.objlist[-1]:
            if self.vb: print ' - Object',self.currentobj,' was the last in the list.\n   Quitting GUI.'
            self.quitting = True
            self.quit_but_cmd()
        else:
            newent = np.where(self.objlist == self.currentobj)[0]+1
            self.currentobj = self.objlist[newent][0]
            self.openpngs()
            self.labelvar.set(self.infostring())

            # load new data for plot and replot
            self.dataPlot_loaddata()
            self.dataPlot_plot(refresh=True,newobj=True)
            self.DPxlow_full, self.DPxhigh_full, self.DPylow_full, self.DPyhigh_full = \
                self.dataPlot_getwindowinfo() # store full window

            if self.fitsauto: # loading fits files automatically
                if self.xpa:
                    self.openfits_but_cmd_xpa()
                else:
                    self.openfits_but_cmd()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def next_but(self,position):
        self.next = Button(self)
        self.next["text"] = "(8) Next object (save)"
        self.next["command"] = self.next_but_cmd
        self.next.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def next_but_cmd(self):
        """
        Command for next button
        """
        self.reset()

        if self.currentobj == self.objlist[-1]:
            if self.vb: print ' - Object',self.currentobj,' was the last in the list.\n   Quitting GUI.'
            self.quitting = True
            self.quit_but_cmd()
        else:
            newent = np.where(self.objlist == self.currentobj)[0]+1
            self.currentobj = self.objlist[newent][0]
            self.openpngs()
            self.labelvar.set(self.infostring())

            # load new data for plot and replot
            self.dataPlot_loaddata()
            self.dataPlot_plot(refresh=True,newobj=True)
            self.DPxlow_full, self.DPxhigh_full, self.DPylow_full, self.DPyhigh_full = \
                self.dataPlot_getwindowinfo() # store full window

            if self.fitsauto: # loading fits files automatically
                if self.xpa:
                    self.openfits_but_cmd_xpa()
                else:
                    self.openfits_but_cmd()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def prev_but(self,position):
        self.prev= Button(self)
        self.prev["text"] = "(7) Previous object"
        self.prev["command"] = self.prev_but_cmd
        self.prev.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def prev_but_cmd(self):
        """
        Command for previous button
        """
        self.reset()

        if self.currentobj == self.objlist[0]:
            if self.vb: print ' - At first object of list...'
        else:
            newent = np.where(self.objlist == self.currentobj)[0]-1
            self.currentobj = self.objlist[newent][0]
            self.openpngs()
            self.labelvar.set(self.infostring())

            # load new data for plot and replot
            self.dataPlot_loaddata()
            self.dataPlot_plot(refresh=True,newobj=True)
            self.DPxlow_full, self.DPxhigh_full, self.DPylow_full, self.DPyhigh_full = \
                self.dataPlot_getwindowinfo() # store full window

            if self.fitsauto: # loading fits files automatically
                if self.xpa:
                    self.openfits_but_cmd_xpa()
                else:
                    self.openfits_but_cmd()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def quit_but(self,position):
        """
        Set up the quit button
        """
        self.QUIT = Button(self)
        self.QUIT["text"] = "QUIT MiG1D"
        self.QUIT["command"] = self.quit_but_cmd
        self.QUIT.grid(row=position[0],column=position[1],columnspan=position[2],sticky=W)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def quit_but_cmd(self):
        """
        Command for quit button
        """
        if self.quitting == False: self.reset() # Only reset if quit_but_cmd was activated by quit button
        self.quitting = True
        self.fout.close()
        self.closewindows()
        self.dataPlotManager.destroy()
        if self.outcheck: self.checkoutput()
        self.quit()
        if self.vb: print ' - Quit MiG1D successfully'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def checkoutput(self):
        """
        Checking the output to see if it is as expected
        """
        data      = np.genfromtxt(self.outfile,comments='#',skip_header=2,names=True)
        Nobjout   = len(np.unique(data['ID']))

        if self.vb: print ' - OUTPUTCHECK: Found '+str(Nobjout)+' objects in output. '+\
                          'Input objlist contained '+str(len(self.objlist))+' objects'
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def keyboard_cmd(self,event):
        """
        Commands for keyboard shortcuts
        """
        cmd = event.char

        focuson = self.focus_get() # check where the focus is
        if  (focuson == self.comments) or (focuson == self.byhandz):
            pass
        else:
            #if self.vb: print '   Keyboard shortcut: ',cmd
            keycmd    = []
            keynames  = []
            keynumber = []
            for ii, key in enumerate(self.keys):
                keycmd.append(key[1])
                keynames.append(key)
                keynumber.append(ii)

            if cmd in keycmd:
                thiskey = keynames[np.where(np.asarray(cmd) == np.asarray(keycmd))[0]]
                if cmd in self.sliders:
                    sliderval = int(self.keys[thiskey].get())
                    if sliderval == 4:
                        self.sliderdic[thiskey].set(0)
                    else:
                        self.sliderdic[thiskey].set(sliderval+1)
                elif cmd in self.empty:
                    pass
                else:
                    self.cbdic[thiskey].toggle()

            elif cmd == 'x':
                self.comments.focus_set()

            elif cmd == 'm':
                sliderval = float(self.slidersmooth.get())
                self.slidersmooth.set(sliderval+0.1)

            elif cmd == 'n':
                sliderval = float(self.sliderz.get())
                self.sliderz.set(sliderval+0.1)

            elif cmd == 'o':
                self.skybox.toggle()

            elif cmd == 'p':
                self.err1Dbox.toggle()

            elif cmd == 'q':
                sliderval = int(self.sliderzqual.get())
                if sliderval == 4:
                    self.sliderzqual.set(0)
                else:
                    self.sliderzqual.set(sliderval+1)

            elif cmd == 'r':
                self.dataPlot_redrawbutton_cmd()

            elif cmd == 's':
                self.dataPlot_savebutton_cmd()

            elif cmd == 't':
                self.mzsbox.toggle()

            elif cmd == 'u':
                self.byhandz.focus_set()

            elif cmd == 'z':
                self.dataPlot_fullzoombutton_cmd()

            elif cmd == '0':
                if self.xpa:
                    self.openfits_but_cmd_xpa()
                else:
                    self.openfits_but_cmd()

            elif cmd == '7':
                self.prev_but_cmd()

            elif cmd == '8':
                self.next_but_cmd()

            else:
                #if self.vb: print ' - Invalid keyboard shortcut :',cmd # prints shift+tab key as well
                pass

#-------------------------------------------------------------------------------------------------------------
#                                                      END
#-------------------------------------------------------------------------------------------------------------

