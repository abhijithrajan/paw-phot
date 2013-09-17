#! /usr/bin/env python 

#=======================================#
#	     	PAW-Phot					#
#    Written by: Paul Anthony Wilson	#
# 	   paw@astro.ex.ac.uk				#
#=======================================#


import sys
from os import getcwd, listdir
import app_phot as ap
import param as par

def phot(folder):
  workdir = getcwd()
  datadir=par.datadir[0]
  dir_contents = listdir(workdir+'/'+folder+'/'+datadir)
  ds9s = [fn for fn in dir_contents if fn.startswith('ds9') and fn.endswith('v2.reg')]
  ds9_name = ds9s[0]
  
  print "\n Loading the photometry packages from IRAF"
  ap.load_phot()
  if par.datadir[0] == "data_binned":
    ap.cd(workdir+'/'+folder+'/'+par.datadir[0])
    ap.app_phot(workdir+'/'+folder,ds9_name)
    ap.cd(workdir)
  else:
    ap.cd(workdir+'/'+folder+'/data')
    ap.app_phot(workdir+'/'+folder)
    ap.cd(workdir)
