#! /usr/bin/env python 

#=======================================#
#    Written by: Paul Anthony Wilson	#
# 	   	  paw@astro.ex.ac.uk	#
#=======================================#

from os  import chdir, path, system, listdir
import os.path
from pyraf import iraf
from pylab import *
import pyraf.iraf as iraf
import sys

import param as par

def fwhm(folder):
  iraf.noao(_doprint=0)
  iraf.obsutil(_doprint=0)
  iraf.unlearn("psfmeasure")
  
  datadir=par.datadir[0]

  name = folder[15:21]
  
  workdir = os.getcwd()
  dir_contents = listdir(workdir+'/'+folder+'/'+datadir)
  ds9s = [fn for fn in dir_contents if fn.startswith('ds9') and fn.endswith('.reg')]
  ds9 = loadtxt(workdir+'/'+folder+'/'+datadir+'/'+ds9s[0])
  
  
  '''
  coords = [fn for fn in dir_contents if fn.endswith('.ctr.1')]
  #coords = [fn for fn in dir_contents if fn.endswith('.coords')]
  coords = sorted(coords)
  if not coords:
    if not os.path.isfile(workdir+'/'+folder+'/data/ds9.reg'):
      print "\nSelect your reference stars first. This is done by creating a ds9.reg file with star coordinates."
      sys.exit()
  '''    
  print "\nMeasuring the fwhm of selected reference stars..."
  iraf.noao.obsutil.psfmeasure.coords = "mark1"
  iraf.noao.obsutil.psfmeasure.display = "no"
  iraf.noao.obsutil.psfmeasure.size = "FWHM"
  iraf.noao.obsutil.psfmeasure.scale = "1"
  iraf.noao.obsutil.psfmeasure.radius = "11"
  iraf.noao.obsutil.psfmeasure.swidth = "15"
  iraf.noao.obsutil.psfmeasure.iterati = "1"
  iraf.noao.obsutil.psfmeasure.imagecu = ds9
  iraf.noao.obsutil.psfmeasure.graphcur = "src/q.reg"
  iraf.noao.obsutil.psfmeasure.logfile = "fwhm.log"

  if os.path.isfile('fwhm.log'):
    os.system('rm fwhm.log')
    
  imgs = [fn for fn in dir_contents if fn.startswith('ali') and fn.endswith('.fits')]
  imgs = sorted(imgs)
 

  for i in range(len(imgs)):
    iraf.noao.obsutil.psfmeasure.imagecu = workdir+'/'+folder+'/'+datadir+'/'+ds9s[0]
    iraf.noao.obsutil.psfmeasure(workdir+'/'+folder+'/'+datadir+'/'+imgs[i])
  
  N_stars = len(ds9)
  fwhm = [[] for _ in range(len(imgs))]
  values = [ line for line in open('fwhm.log') if '.' in line and 'NOAO' not in line]
  j = 0
  for i in range(len(values)):
    if values[i][2:9] == 'Average':
      j += 1
    if values[i][2:9] != 'Average': 
      fwhm[j].append(float(values[i][41:47]))

  if datadir=="data_binned":
    f = open('fwhm_binned.txt', 'w+')
  else:
    f = open('fwhm.txt', 'w+')
  for i in range(len(imgs)):
    fwhm[i] = array(fwhm[i])
    print >> f, np.median(fwhm[i]), imgs[i]
  f.close()

  
  if datadir=="data_binned":
    os.system('mv fwhm_binned.txt '+folder+'/'+datadir)
    os.system('mv fwhm.log '+folder+'/'+datadir)
  else:
    os.system('mv fwhm.txt '+folder+'/'+datadir)
    os.system('mv fwhm.log '+folder+'/'+datadir)
