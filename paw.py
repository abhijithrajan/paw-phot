#! /usr/bin/env python 

#=======================================#
#  PAW-Phot v1.0 Beta
#  Written by: Paul Anthony Wilson
#  University of Exeter
#  paw@astro.ex.ac.uk
#=======================================#

import os
from os import getcwd, chdir, path, makedirs
import sys
import numpy as np
from operator import itemgetter

import pyfits
from pyfits import getheader

from src import unix as u
from src import param as par
from src import fwhm,phot,plot,best_aperture

# Start Program
print '\n PAW-Phot'
print ' Version: 1.0 Beta'
print ' paw@astro.ex.ac.uk'

def show_commands():
  print "\n\t\t"+u.bold('PAW-Phot Commands:')+"\n"
  print "\t"+u.green('fwhm')+"\tMeasures the fwhm of selected stars in field."
  print "\t\ti.e. "+u.green('fwhm 01.')#+" or "+u.green('fwhm all')+"."
  print "\t"+u.green('phot')+"\tPerforms the photometry."
  print "\t"+u.green('plot')+"\tCreates the light curves."
  print "\t"+u.red(' q')+"\tQuit the program.\n"

show_commands()

# Create a login.cl file if it does not already exist
if not os.path.isfile('login.cl'):
  os.system("mkiraf -t xgterm -i -q")
  print "\n No login.cl and uparm directory found. It has now been created.\n"

input=1
while 1 :
    input = raw_input("paw >> ")
    num = input[5:7]
    rootdir = par.rootdir[0]
    datadir=rootdir+'/'+par.datadir[0]
    dir_contents = os.listdir(datadir)
    folders = [s for s in dir_contents if s.startswith("OBJ")]
    folders = sorted(folders)

    if input in ['exit','q','quit']:
        exit()
    
    # FWHM
    # ============================================= #
    elif 'fwhm' in input:
        if input[5:8] == 'all':
          print "Doing all"
          for num in range(len(folders)):
            print "\nWORKING WITH: ",folders[int(num)],"\n"
            fwhm.fwhm(datadir+'/'+folders[int(num)],folders)
	else:
	  	  fwhm.fwhm(datadir+'/'+folders[int(num)-1],folders)
    # ============================================= #


    # PHOTOMETRY
    # ============================================= #
    elif 'phot' in input:
        if input[5:8] == 'all':
          print "Doing all"
          for num in range(len(folders)):
            phot.phot(datadir+'/'+folders[int(num)],folders)
        else:       
          phot.phot(datadir+'/'+folders[int(num)-1],folders)
    # ============================================= #
    

    # PLOT and LIGHT CURVE CREATION
    # ============================================= #
    elif 'plot' in input:
        f = open('data.txt', 'w+')
        f.write("# OBJ\t\tChi2\tChi2_m\tRobust\tSTD_target\tSTD_Refs:\tSTD_Master:\tInital\tBad_ref\tRefs\tDOF\tA\tA_err\tp-val\tPoss_var\n")
        f.close()
        if input[5:8] == 'all':
          print "Doing all"
          dir_contents = os.listdir('all_data')
          folders = [s for s in dir_contents if s.startswith("OBJ")]
          folders = sorted(folders)
          for num in range(len(folders)):
			try:
			  if num < 9 and num+1:# not in [6]:
				print "\nWORKING WITH: ",datadir+'/'+folders[int(num)-1],folders,"\n"
				plot.plot(datadir+'/'+folders[int(num)-1],folders)
			  else:
				#sys.exit()
				if num+1 not in par.ignore_objects and num+1 > -45:
				  print "\nWORKING WITH: ",datadir+'/'+folders[int(num)-1],folders,"\n"
				  plot.plot(datadir+'/'+folders[int(num)-1],folders)
			except IndexError:
			  print u.yellow("Trying again...(IndexError)")
			  if num < 9 and num+1:# not in [6]:
				print "\nWORKING WITH: ",datadir+'/'+folders[int(num)-1],folders,"\n"
				plot.plot(datadir+'/'+folders[int(num)-1],folders)
			  else:
				#sys.exit()
				if num+1 not in par.ignore_objects and num+1 > -45:
				  print "\nWORKING WITH: ",datadir+'/'+folders[int(num)-1],folders,"\n"
				  plot.plot(datadir+'/'+folders[int(num)-1],folders)
        else:
	  	  plot.plot(datadir+'/'+folders[int(num)-1])
    # ============================================= #
    
    elif 'help' in input:
        show_commands()
    else:
        print "Command not recognized"
