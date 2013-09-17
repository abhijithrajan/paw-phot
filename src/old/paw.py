#! /usr/bin/env python 

#=======================================#
#	PAW-Phot v2.0 Beta		#
#  Written by: Paul Anthony Wilson	#
#	paw@astro.ex.ac.uk		#
#=======================================#

import os
from os import getcwd, chdir, path, makedirs
import time, sys
import numpy as np
from operator import itemgetter

import pyfits
from pyfits import getheader				# Jobb med

from src import unix as u
from src import fwhm
from src import phot
from src import plot
from src import best_aperture

# Start Program
print '\n PAW-Phot'
print ' Version: 2.0 Beta'
starttime=time.ctime()
print ' paw@astro.ex.ac.uk'

def show_commands():
  print "\n\t\t"+u.bold('PAW-Phot Commands:')+"\n"
  print "\t"+u.green('fwhm')+"\tMeasures the fwhm of selected stars in field."
  print "\t\ti.e. "+u.green('fwhm 27')+" or "+u.green('fwhm all')+"."
  print "\t"+u.green('phot')+"\tPerforms the photometry."
  print "\t"+u.green('plot')+"\tCreates the light curves."
  print "\t"+u.green('bestapp')+"\tFinds the optimal aperture."
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
    dir_contents = os.listdir('all_data')
    folders = [s for s in dir_contents if s.startswith("OBJ")]
    folders = sorted(folders)
    # Python 3 users
    # input = input(">> ")
    if input in ['exit','q','quit']:
        exit()
    # FWHM
    elif 'fwhm' in input:
        if input[5:8] == 'all':
	  print "Doing all"
	  dir_contents = os.listdir('all_data')
	  folders = [s for s in dir_contents if s.startswith("OBJ")]
	  folders = sorted(folders)
	  for i in range(len(folders)):
	    print "Performing FWHM measurements for all data sets\t"+str(round(i/float(len(folders))*100,2))+"% done"
            if not os.path.isfile('all_data/OBJ'+str(i+1)+'/data/fwhm.txt'):
	      print "Working in: ",'all_data/OBJ'+str(i+1)
	      fwhm.fwhm('all_data/OBJ'+str(i+1))
	else:
	  	  fwhm.fwhm('all_data/OBJ'+num+folders[int(num)-1][5:])

    # PHOTOMETRY
    elif 'phot' in input:
        if input[5:8] == 'all':
	  print "Doing all"
	  dir_contents = os.listdir('all_data')
	  folders = [s for s in dir_contents if s.startswith("OBJ")]
	  folders = sorted(folders)
	  # del folders[18]
	  for i in range(len(folders)):
	    print "Performing FWHM measurements for all data sets\t"+str(round(i/float(len(folders))*100,2))+"% done"
            if not os.path.isfile('all_data/OBJ'+str(i+2)+'/data/object.list'):
	      print "Working in: ",'all_data/OBJ'+str(i+2)
	      phot.phot('all_data/OBJ'+str(i+2))
	else:
          phot.phot('all_data/OBJ'+num+folders[int(num)	-1][5:])

    # PLOT
    elif 'plot' in input:
        if input[5:8] == 'all':
          print "Doing all"
          dir_contents = os.listdir('all_data')
          folders = [s for s in dir_contents if s.startswith("OBJ")]
          folders = sorted(folders)
          for i in range(len(folders)):
			''' Dirty header code which can be removed
			print i,folders[i],'OBJ'+str(i+1)
		
			#print i,folders[i],'all_data/OBJ'+str(i+1)
			img_contents = os.listdir('all_data/'+folders[i]+'/data/')
			images = [s for s in img_contents if s.endswith(".fits")]
			headers = []
			headers.append(pyfits.open('all_data/'+folders[i]+'/data/'+images[0]))
			print folders[i],headers[0][0].header['HIERARCH ESO INS FILT1 NAME']
			'''

			#if not os.path.isfile('all_data/OBJ'+str(i+1)+'/ref_stars.list'):
			#print 'all_data/OBJ'+str(i+1)
			#if (i+1) not in [0,14,17,26,27,31,37,41,44,45,46,57,60,62,63,69]:
			if (i+1) not in [26,27]:
			  plot.plot('all_data/OBJ'+str(i+1))#+input[5:7])
        else:
	  	  plot.plot('all_data/OBJ'+num+folders[int(num)	-1][5:])
		  '''
			  for i in range(70):
			if i not in [0,14,17,26,27,31,37,41,44,45,46,57,60,62,63,69]:#in range(69):    #Check 27    
				  #plot.plot('all_data/OBJ'+input[5:7]) #45 variable? 57? 63 is variable I think
				  #if i > 67:
			  plot.plot('all_data/OBJ'+str(i))
		  '''

    # BESTAPP
    elif 'bestapp' in input:
    	if input[8:11] == 'all':
    	  print "Doing all"
    	  dir_contents = os.listdir('all_data')
          folders = [s for s in dir_contents if s.startswith("OBJ")]
          folders = sorted(folders)
          for i in range(len(folders)):
            best_aperture.best_app('all_data/'+folders[i]+'/not_used','off')
    	else:
          num = input[8:10]	# Change due to bestapp word being longer than fwhm,phot,plot
          best_aperture.best_app('all_data/OBJ'+num+folders[int(num)	-1][5:],'on')        
    elif 'help' in input:
        show_commands()
    else:
        print "Command not recognized"


