#! /usr/bin/env python 

#=======================================#
#	     	PAW-Phot					#
#    Written by: Paul Anthony Wilson	#
# 	   paw@astro.ex.ac.uk				#
#=======================================#


import sys
from os import getcwd, listdir
import a_plot as p
import param as par
from pylab import loadtxt


def plot(folder):
  workdir = getcwd()
  datadir=par.datadir[0]
  dir_contents = listdir(workdir+'/'+folder+'/'+datadir)
  ds9s = [fn for fn in dir_contents if fn.startswith('ds9') and fn.endswith('v2.reg')]
  ds9 = loadtxt(workdir+'/'+folder+'/'+datadir+'/'+ds9s[0])
  print "\n Creating data...\n"
  p.do_plot(workdir+'/'+folder,ds9)
