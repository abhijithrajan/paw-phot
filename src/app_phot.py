#! /usr/bin/env python 

#=======================================#
#    Written by: Paul Anthony Wilson	#
# 	   	  paw@astro.ex.ac.uk			#
#=======================================#

import os, sys
import shutil
from os import getcwd,listdir,chdir
import numpy as np
from pylab import loadtxt

import pyfits
from pyraf import iraf

import param as par

def cd(target):
    chdir(target)
    iraf.cd(target)

def mkdir(dirname):
   "Creates a directory if it didn't exist before"
   if not os.path.isdir(dirname):
      os.mkdir(dirname)

def load_phot():
   iraf.digiphot(_doprint=0)
   iraf.daophot(_doprint=0)
   
def ref_stars():
   "Creates a file containg the star coordinates (see param.py) if coordinate files are not found."
   dir_contents = listdir(workdir+'/'+datadir)
   coords = [fn for fn in dir_contents if fn.endswith('.ctr.1')]
   #coords = [fn for fn in dir_contents if fn.endswith('.coords')]
   coords = sorted(coords)
   if not coords:
     if not os.path.isfile('ds9.reg'):
       print "How many reference stars are you using?"
       N_ref_stars = int(raw_input())
       ds9 = open('ds9.reg', 'w+')
       print >> ds9, par.target_star_loc[0],par.target_star_loc[1]
       for i in range(N_ref_stars):
         print >> ds9, par.ref_x[i],par.ref_y[i]
       ds9.close()

def set_phot_params(num,folder):
   global fwhm
   datadir = par.datadir[0]
   if datadir=='data_binned':
     fwhm = loadtxt(folder+'/'+datadir+'/fwhm_binned.txt',usecols=(0,))
   else:
     fwhm = loadtxt(folder+'/'+datadir+'/fwhm.txt',usecols=(0,))

   # Set centerpars
   iraf.centerpars.setParam('calgori','centroid')
   iraf.centerpars.setParam('cbox',7)
   iraf.centerpars.setParam('maxshift',5)

   # Set datapars
   iraf.datapars.setParam('fwhmpsf',fwhm[num])
   iraf.datapars.setParam('sigma',7.8)
   #iraf.datapars.setParam('datamin',1000)
   #iraf.datapars.setParam('datamax',45000)
   iraf.datapars.setParam('epadu',5.385)
   iraf.datapars.setParam('readnoi',2.1)
   iraf.datapars.setParam('exposure','EXPTIME')
   #iraf.datapars.setParam('airmass','HIERARCH ESO TEL AIRM START')
   #iraf.datapars.setParam('filter','HIERARCH ESO INS FILT1 NAME')
   #iraf.datapars.setParam('obstime','READTIME')

   # Set findpars
   iraf.findpars.setParam('thresho',3) # Threshold in sigma for feature detection
   iraf.findpars.setParam('nsigma',1.5) # Width of convolution kernel in sigma - whatever that means...

   # Set fitskypars
   iraf.fitskypars.setParam('salgorithm','centroid')
   #iraf.fitskypars.setParam('annulus',12) # Inner radius of sky annulus in scale units
   #iraf.fitskypars.setParam('dannulus',15) # Width of sky annulus in scale units

   # Saving parameters in files
   iraf.centerpars.saveParList(filename='center.par')
   iraf.datapars.saveParList(filename='data.par')
   iraf.fitskypars.saveParList(filename='fitsky.par')
   iraf.photpars.saveParList(filename='phot.par')

   # Set some IRAF settings
   iraf.set(imtype="fits")
 
   return fwhm


def app_phot(folder,ds9_name):

   print '____________________________________________________'

   if os.path.isfile("fitsky.par"):
     print "\n Setting initial parameter files		Done"
   else:
     print "\n Please supply the following information:\n"
  
   datadir=par.datadir[0]
   
   # Making a list of image files to be used:
   dir_contents = listdir(folder+'/'+datadir)
   #imgs = [fn for fn in dir_contents if fn.startswith('rfb') and fn.endswith('.fits')]
   imgs = [fn for fn in dir_contents if fn.startswith('ali') and fn.endswith('.fits')]
   #imgs = [fn for fn in dir_contents if fn.startswith('shft') and fn.endswith('.fits')]
   imgs = sorted(imgs)
   #for i in range(len(imgs)):
   #  print imgs[i]
   #sys.exit()
   #coords = [fn for fn in dir_contents if fn.endswith('.ctr.1')]
   #coords = [fn for fn in dir_contents if fn.endswith('.coords')]
   #coordinates = sorted(coords)
   #bad_pix = [fn for fn in dir_contents if fn.startswith('bad_pix_mask')]
   
   print "Creating object list"
   
   print folder+'/'+datadir
   photometry_performed = len([fn for fn in dir_contents if 'xyf.list' in fn]) > 0
   
   if photometry_performed:
   
      print "\n Photometry					Done"
      print '____________________________________________________\n\n'

   else:

      # Determine median FWHM for this dataset
      fwhm = []
      if datadir=="data_binned":
        fwhm2 = open(folder+'/'+datadir+'/fwhm_binned.txt', 'r')
      else:
        fwhm2 = open(folder+'/'+datadir+'/fwhm.txt', 'r')
      eachfwhm = fwhm2.readlines()
      fwhm2.close()
 
      for eachline in eachfwhm:
         p = eachline.split()
         fwhm.append(float(p[0]))

      seeing_value = np.median(fwhm)
      print seeing_value
      
      variable_aperture=par.variable_aperture[0]

      # First aperture radius
      print "\nFirst Aperture Radius in pixels:?"
      if variable_aperture=="no":
        Ap1 = 1.5*seeing_value   #raw_input("(e.g. 5.0): ")
      else:
        Ap1 = 1.5

      # Last aperture radius
      print "Last Aperture Radius in pixels:?"
      if variable_aperture=="no":
        Ap2 = 1.5*seeing_value   #raw_input("(e.g. 5.0): ")
      else:
        Ap2 = 1.5

      # First aperture radius
      print "First Outer Aperture Radius (pixels)?"
      if variable_aperture=="no":
        Ap3 = 1.5*seeing_value   #raw_input("(e.g. 5.0): ")
      else:
        Ap3 = 1.5

      # Last aperture radius
      print "Last Outer Aperture Radius (pixels)?"
      Ap4 = 3.0*seeing_value	#22.  #raw_input("(e.g. 60.0): ")

      # First aperture radius
      print "Starting Width of Outer Aperture Radius (pixels)?"
      Ap5 = 10#raw_input("(e.g. 10.0): ")

      # Last aperture radius
      print "End With of Last Outer Aperture Radius (pixels)?"
      Ap6 = 10#raw_input("(e.g. 15.0): ")

      # step apertures in 0.5 pixel steps
      #NumAps=np.arange(float(Ap1),float(Ap2)+0.1,0.1)
      if Ap1 == Ap2:
        a_range = [Ap1]
        skyrad = [Ap3]
      else:
        a_range = np.arange(float(Ap1),float(Ap2)+0.25,0.25) 	# Step size
        skyrad=np.arange(float(Ap3),float(Ap4)+1.0,1.0)
      anwidth=np.arange(float(Ap5),float(Ap6)+5.0,5.0)

   for j in range(0,len(a_range)):
      for k in range(0,len(skyrad)):
          for l in range(0,len(anwidth)):
              if variable_aperture=="no":
                if datadir=="data_binned":
                  newsubdir = "/B_binned_" + str(a_range[j]) + "xFWHM"# + str(skyrad[k]) + "D" + str(anwidth[l])
                else:
                  newsubdir = "/B" + str(a_range[j]) + "xFWHM_S" + str(skyrad[k]) + "D" + str(anwidth[l])
              else:
                if datadir=="data_binned":
                  newsubdir = "/V_binned_" + str(a_range[j]) + "xFWHM"# + str(skyrad[k]) + "D" + str(anwidth[l])
                else:
                  newsubdir = "/V" + str(a_range[j]) + "xFWHM_S" + str(skyrad[k]) + "D" + str(anwidth[l])
              newdir=str(folder)+str(newsubdir)

              print "REMOVING"            
              os.system('rm -r ../V*')
              
              if os.path.isdir(newdir):
                shutil.rmtree(newdir)
                os.mkdir(newdir)
              else:
                os.mkdir(newdir)
     
              # Set outer photpars
              iraf.fitskypars.setParam('annulus',skyrad[k]) # Inner radius of sky annulus in scale units
              iraf.fitskypars.setParam('dannulus',anwidth[l]) # Width of sky annulus in scale units

              iraf.phot.setParam('interactive','no')
              iraf.phot.setParam('verify','no')
              iraf.phot.setParam('verbose','yes')

     
              # Set inner photpars
              for m in range(len(imgs)):
                set_phot_params(m,folder)
                if Ap1 == Ap2:
                  NumAps=[fwhm[m]*float(Ap1)]
                  skyrad = [fwhm[m]*float(Ap3)]
                else:
                  NumAps=fwhm[m]*a_range
		
                if variable_aperture=="no":
                  iraf.photpars.setParam('apertur',a_range[j]) # List of aperture radii in scale units
                  print '\n Aperture: ' + str(a_range[j]) + '\tOuter radius: ' + str(skyrad[k]) + '\tThickness: '+str(anwidth[l])+'\n'
                else:
                  iraf.photpars.setParam('apertur',NumAps[j]) # List of aperture radii in scale units
                  iraf.fitskypars.setParam('annulus',skyrad[k]) # Inner radius of sky annulus in scale units
                  print '\n Aperture: ' + str(NumAps[j]) + '\tOuter radius: ' + str(skyrad[k]) + '\tThickness: '+str(anwidth[l])+'\n'
                
                print "Doing: ",j+1,"/",len(a_range),"\t\t\t",k+1,"/",len(skyrad),"\t\t\t",l+1,"/",len(anwidth)
                print round((float(m)/len(imgs))*100.,1),"% done"
                iraf.photpars.saveParList(filename="phot-loop.par")
                iraf.phot.saveParList(filename="photpars.pars")

                iraf.daophot.phot(

                ####	OBJ parameters	     ####
                str(folder+'/'+datadir+'/'+imgs[m]),
                coords=folder+'/'+datadir+"/"+ds9_name,
                #coords=str(coordinates[m]),
                output="default",
                datapar="data.par",
                centerp="center.par",
                fitskyp="fitsky.par",
                photpar="phot-loop.par",
                interac="no",
                verify="no",
                update="no",
                verbose="no")
                #################################

                
                print '___________________________________________________________________\n'      
                # Dumping data from the .mag.1 files into text files
              iraf.daophot.pdump(	# Getting the position and flux of stars

              ####	OBJ parameters	     ####
              					#
              "*.mag.1",			#
              fields="XCEN,YCEN,FLUX,MAG,MERR,MSKY,XCENTER,YCENTER",#
              expr="yes",			#
              headers="no",			#
              Stdout=folder+'/'+datadir+'/'+'xyf.list')		#
              					#
              #################################

              # Creating directories
              mag = newdir+"/"+"mag"
              mkdir(mag)
              param = newdir+"/"+"param"
              mkdir(param)                  
                  
              os.system( 'mv *.mag.1 ' +mag  )
              os.system( 'mv '+folder+'/'+datadir+'/'+'xyf.list'+' '+newdir  )
              os.system( 'mv *par* ' +param  )
            
   print "Median seeing:\t",seeing_value
   print "Photometry done\n"
              
   return
