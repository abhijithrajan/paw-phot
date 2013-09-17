#! /usr/bin/env python 

#=======================================#
#    Written by: Paul Anthony Wilson	#
# 	  paw@astro.ex.ac.uk		#
#=======================================#

import sys
import pyfits
from pyfits import getheader
from numpy import loadtxt, array
import numpy as np
from pylab import *
from os import listdir#, path
import os
import time
import matplotlib
import math
import random
from scipy import stats

import unix as u
import param as par
import bin_calc as bc
import calc as c

ion()

def get_airmass_JD(folder):
  global JD, imgs
  datadir=par.datadir[0]
  dir_contents = listdir(folder+'/'+datadir)
  imgs = [s for s in dir_contents if s.startswith("ali") and s.endswith(".fits")]
  imgs = sorted(imgs)
  headers = []
  airmass = []
  JD = []
  for i in range(len(imgs)):
    headers.append(pyfits.open(folder+'/'+datadir+'/'+imgs[i]))
    JD.append(headers[i][0].header['MJD-COMB'])
    airmass.append(headers[i][0].header['HIERARCH ESO TEL AIRM START'])
  return airmass, JD

def do_plot(folder,ds9):

	datadir = par.datadir[0]
	
	''' Defining the number of reference stars '''
	#ds9 = loadtxt(folder+'/'+datadir+'/ds9_'+folder[33:39]+'.reg')
	
	if datadir == "data_binned":
	  fwhm = loadtxt(folder+'/'+datadir+'/fwhm_binned.txt',usecols=(0,))
	else:
	  fwhm = loadtxt(folder+'/'+datadir+'/fwhm.txt',usecols=(0,))
	N_ref_stars = len(ds9)-1

	bad_mag_refs = []
	mag_refs = []	 

	MAD_factor = par.cut_offs[0]	#cut off in the bins

	''' Read Julian Date and airmass info from .fits headers '''
	airmass, JD = get_airmass_JD(folder)
	airmass = np.array(airmass)
	JD = np.array(JD)
	turbulence_fwhm = np.ones(len(fwhm)) * np.array([random.random() for _ in xrange(len(fwhm))])*1e-12
	turbulence_airmass = np.ones(len(airmass)) * np.array([random.random() for _ in xrange(len(airmass))])*1e-12
	fwhm = fwhm+turbulence_fwhm		# A tiny pertubation is added as some of the FWHM values measured by IRAF are identical causing more FWHM elements to be revomved than what is desired.
	airmass = airmass+turbulence_airmass

	''' Calculate the seeing '''
	pixel_scale = par.pixel_scale[0]
	seeing = fwhm*pixel_scale

	''' Establishing directories '''
	dir_contents = listdir(folder)
	if datadir=="data_binned":
	  if par.variable_aperture[0]=="yes":
	    folders = [s for s in dir_contents if s.startswith("V_binned_")]	  
	  else:
	    folders = [s for s in dir_contents if s.startswith("B_binned_")]
	  data_folder = [s for s in dir_contents if s.startswith("data_binned")]
	else:
	  if par.variable_aperture[0]=="yes":
	    folders = [s for s in dir_contents if s.startswith("V")]	  
	  else:
	    folders = [s for s in dir_contents if s.startswith("B")]
	data_folder = [s for s in dir_contents if s.startswith("da") and s.endswith("ta")]
	directory = folder+'/'+folders[0]
	data_directory = folder+'/'+data_folder[0]
	T = loadtxt(directory+'/xyf.list')

	''' Convert JD days --> hours '''
	JD_1 = JD[0]
	JD = (JD-JD_1)#(JD-JD_1)*24

	''' Creating 2D lists list[star][frame] '''
	flux_r = [[] for _ in range(N_ref_stars)]
	mag_r = [[] for _ in range(N_ref_stars)]
	err_r = [[] for _ in range(N_ref_stars)]
	err_r_MAG = [[] for _ in range(N_ref_stars)]
	flux_r_norm = [[] for _ in range(N_ref_stars)]
	err_r_norm = [[] for _ in range(N_ref_stars)]
	stddev = [[] for _ in range(N_ref_stars)]
	stddev_norm = [[] for _ in range(N_ref_stars)]
	variance_norm = [[] for _ in range(N_ref_stars)]

	''' Obtaining the photometry data '''
	x, y, flux, mag, merr, msky, xcenter, ycenter = loadtxt(directory+'/xyf.list', unpack=True)
	x_target = x[::N_ref_stars+1]
	y_target = y[::N_ref_stars+1]

	''' Reading the Target Info '''
	flux_T = flux[::N_ref_stars+1]
	flux_T_MAG = mag[::N_ref_stars+1]

	''' Normalising Target Flux '''
	flux_target_norm = flux_T/np.median(flux_T)#flux_T.mean()
	err_T = merr[::N_ref_stars+1]*flux_T/1.0857	# Convert from mag error to flux error
	err_T_MAG = merr[::N_ref_stars+1]+1e-6
	err_target_norm = err_T/np.median(flux_T)#flux_T.mean()

	iteration = 0
	N_bad_refs = 0

	norm_ref_star_flux_weight_sum_mean_Rs = []
	REF_LCs = []
	err_ref_Rs = []
	mag_refs = []
	''' Reading reference star info and normalising '''
	for i in range(N_ref_stars):
	  flux_r[i] = np.array(flux[i+1::N_ref_stars+1])
	  mag_r[i] = np.array(mag[i+1::N_ref_stars+1])
	  err_r[i] = np.array(merr[i+1::N_ref_stars+1]*flux_r[i]/1.0857)
	  err_r_MAG[i] = np.array(merr[i+1::N_ref_stars+1] )
	  ''' Normalising the reference flux and errors '''
	  flux_r_norm[i] = np.array(flux_r[i]/np.median(flux_r[i]))
	  # Adding 1e-12 as IRAF has errors stopping at 0.0001 and I don't want to divide by zero:
	  err_r_norm[i] = np.array(err_r[i]/np.median(flux_r[i])+1e-12)#/flux_r[i].mean()+1e-12)
	  stddev[i] = flux_r[i].std()
	  stddev_norm[i] = flux_r_norm[i].std()
	  variance_norm[i] = array(err_r_norm[i]*err_r_norm[i])
	      
	  mag_refs.append(np.median(mag_r[i]))


	while (iteration < 2):

	  N_bad_refs = 0
	  print "\n\nITERATION: ",iteration
	  print 'Ref stars: ',N_ref_stars
	  bad_ref_list = []
	  if (iteration == 0 and N_ref_stars > par.cut_offs[2]):
 	    diff = []
 	    sim_refs = []
 	    for i in range(N_ref_stars):
 	      diff.append( [abs(np.median(mag_refs[i])-np.median(flux_T_MAG)),i] )
 	    diff = sorted(diff)
 	    for i in range(par.cut_offs[2]):
 	      sim_refs.append(diff[i][1])
 	    dissimilar_refs = [x for x in np.arange(0,N_ref_stars) if x not in sim_refs]	    
	    print 'Removing ',len(dissimilar_refs),' refs'	    
	    bad_ref_list = dissimilar_refs

	    # Deleting bad refs
	    flux_r = [flux_r[i] for i in sim_refs]
	    mag_r = [mag_r[i] for i in sim_refs]
	    err_r = [err_r[i] for i in sim_refs]
	    err_r_MAG = [err_r_MAG[i] for i in sim_refs]
	    flux_r_norm = [flux_r_norm[i] for i in sim_refs]
	    err_r_norm = [err_r_norm[i] for i in sim_refs]
	    stddev = [stddev[i] for i in sim_refs]
	    stddev_norm = [stddev_norm[i] for i in sim_refs]
	    variance_norm = [variance_norm[i] for i in sim_refs]


	  
	  elif (os.path.isfile('temp.txt') and len(np.loadtxt('temp.txt',usecols=(0,),ndmin=1))>0. and N_ref_stars>par.cut_offs[3] and iteration == 1):
	    if (iteration == 1 and N_ref_stars > par.cut_offs[3]):
 	      diff = []
 	      sim_refs = []
 	      for i in range(N_ref_stars):
 	        diff.append( [abs(np.median(mag_refs[i])-np.median(flux_T_MAG)),i] )
 	      diff = sorted(diff)
 	      for i in range(par.cut_offs[2]):
 	        sim_refs.append(diff[i][1])
 	      dissimilar_refs = [x for x in np.arange(0,N_ref_stars) if x not in sim_refs]	    
 	      print 'Removing ',len(dissimilar_refs),' refs'	    
 	      bad_ref_list = dissimilar_refs
	    

	    bad_ref_stars = np.loadtxt('temp.txt',usecols=(0,),ndmin=1)

	    N_bad_refs = len(bad_ref_stars)
	    
	    bad_ref_list = []
	    print 'Removing ',len(bad_ref_stars),' refs based on MAD'
	    for i in range(len(bad_ref_stars)):
	      bad_ref_list.append(int(bad_ref_stars[i]))

	    best_refs = [x for x in np.arange(0,N_ref_stars) if x not in bad_ref_list]

	    # Deleting bad refs based on MAD
	    flux_r = [flux_r[i] for i in best_refs]
	    mag_r = [mag_r[i] for i in best_refs]
	    err_r = [err_r[i] for i in best_refs]
	    err_r_MAG = [err_r_MAG[i] for i in best_refs]
	    flux_r_norm = [flux_r_norm[i] for i in best_refs]
	    err_r_norm = [err_r_norm[i] for i in best_refs]
	    stddev = [stddev[i] for i in best_refs]
	    stddev_norm = [stddev_norm[i] for i in best_refs]
	    variance_norm = [variance_norm[i] for i in best_refs]
   
	    
	  N_ref_stars = N_ref_stars - len(bad_ref_list)      

	  
	  ''' Find the mean of the raw flux of the reference stars '''
	  norm_ref_star_flux_sum = np.array([sum(a) for a in zip(*(flux_r_norm))])
	  norm_ref_star_flux_sum = norm_ref_star_flux_sum/N_ref_stars


	  ''' Initialising the 2D arrays to be used in the weighted mean calculations'''
	  weights = [[] for _ in range(N_ref_stars)]
	  weighted_mean = [[] for _ in range(N_ref_stars)]
	  norm_ref_star_flux_weight_sum_mean = [[] for _ in range(N_ref_stars)]

	  ''' Calculating the weighted mean '''
	  for i in range(N_ref_stars):
	    weights[i] = 1. / variance_norm[i]		# variance_norm = err_r_norm^2
	    weighted_mean[i] = weights[i]*flux_r_norm[i]
	    
	  weights_sum = np.array([mean(a) for a in zip(*(weights))])
	  weighted_mean_sum = np.array([mean(a) for a in zip(*(weighted_mean))])
	  weighted_mean_sum_norm = weighted_mean_sum/weighted_mean_sum.sum()	# Does sum to 1
	  ''' Creating Master Ref Curve with associated errors '''
	  norm_ref_star_flux_weight_sum_mean = weighted_mean_sum / weights_sum
	  sigma_weighted_mean_ref_stars = np.sqrt(1./(weights_sum*N_ref_stars))

	  ''' Calculating the weighted mean excluding the ref star used
		  to create ref star light curves. i.e.
		  REF_LCs[2] = flux_r_norm[2]/norm_ref_star_flux_weight_sum_mean_Rs[2]
		  where norm_ref_star_flux_weight_sum_mean_Rs[2] does not include the
		  ref star with index 2.
		  '''	  
	  for i in range(N_ref_stars):
	    weights_Rs = weights
	    weighted_mean_Rs = weighted_mean
	      
	    weights_Rs = [x for x in weights_Rs if x not in weights_Rs[i]] # Removing the contribution by the comparison star
	    weights_sum_Rs = np.array([mean(a) for a in zip(*(weights_Rs))])
	    weighted_mean_Rs = [x for x in weighted_mean_Rs if x not in weighted_mean_Rs[i]]
	    weighted_mean_sum_Rs = np.array([mean(a) for a in zip(*(weighted_mean_Rs))])
	      
	    weighted_mean_sum_norm_Rs = weighted_mean_sum_Rs/weighted_mean_sum_Rs.sum()

	    norm_ref_star_flux_weight_sum_mean_Rs.append(weighted_mean_sum_Rs / weights_sum_Rs)	  
	    REF_LCs.append(flux_r_norm[i]/norm_ref_star_flux_weight_sum_mean_Rs[i])
	    err_ref_Rs.append(np.sqrt(err_r_norm[i]**2 + np.median(sigma_weighted_mean_ref_stars)**2))	# IRAF errors

	  master_ref = np.array([np.median(a) for a in zip(*(REF_LCs))])
	  master_ref = master_ref/np.median(master_ref)

	  ''' FINAL ERROR BARS using IRAF values '''
	  err_target_w_mean_norm = np.sqrt(err_target_norm**2 + np.median(sigma_weighted_mean_ref_stars)**2) # Target errors and error of combined ref LC added in quadrature.

	  ''' Sum of reference star errors (as given by IRAF)  '''
	  err_r_norm_sum = np.array([np.median(a) for a in zip(*(err_r_norm))])/np.sqrt(N_ref_stars)		# Errors not using a weigthed mean.

	  ''' Check if in fact, the weighted mean as improved the STDDECV of LC. '''
	  if (iteration == 0):	
		if 	((flux_target_norm/norm_ref_star_flux_weight_sum_mean).std() <= (flux_target_norm/norm_ref_star_flux_sum).std()):
			print "The weighted mean has ",u.green("improved")," the STDDEV of the LC.\n"
		else:
			print "Weighted mean ",u.red("not")," really helping here.\n"
			print (flux_target_norm/norm_ref_star_flux_weight_sum_mean).std(),"!<",(flux_target_norm/norm_ref_star_flux_sum).std(),"\n"
	
		if (sigma_weighted_mean_ref_stars.sum() <= err_r_norm_sum.sum() ):
			print "The weighted mean has ",u.green("improved")," the errors.\n"
		else:
			print "Weighted mean ",u.red("not")," really helping the errors here.\n"

	  ''' The non detrended target light curve '''
	  T_LC = flux_target_norm/norm_ref_star_flux_weight_sum_mean
	  T_LC = T_LC/np.median(T_LC)
	  
	  ''' The non detrended but weighted sample reference star light curve '''
	  #REF_LC = REF_LCs[ref_number]
	
	  if (iteration == 0):
	    # Detrending
	    
	    # AIRMASS
	    A,B,C = np.polyfit(airmass,master_ref, 2) # Fitting a 2D polynomial to the airmass
	    np_airmass = np.arange(airmass.min(),airmass.max(),1e-3)
	    poly2d_AIRMASS = A*np_airmass**2+B*np_airmass+C
	    airmass_effect = A*airmass**2+B*airmass+C
	 
	  if (iteration == 1 and par.variable_aperture[0]=="no"):   
	    # FWHM
	    slope_FWHM, intercept_FWHM, r_value_AIRMASS, p_value_AIRMASS, std_err_AIRMASS = stats.linregress(fwhm,master_ref)
	    print "AIRMASS params:\t",round(A,5),round(B,5),round(C,5)
	    print "FWHM params:\t",round(slope_FWHM,5),round(intercept_FWHM,5)
	    fwhm_effect = slope_FWHM*fwhm+intercept_FWHM
	
	  if (iteration == 1):
	    T_LC = T_LC/airmass_effect
	    if (par.variable_aperture[0]=="no"):
	      T_LC = T_LC/fwhm_effect
	    T_LC = T_LC/np.median(T_LC)
	    
	    for i in range(N_ref_stars):
	      REF_LCs[i] = REF_LCs[i]/airmass_effect
	      if (par.variable_aperture[0]=="no"):
	        REF_LCs[i] = REF_LCs[i]/fwhm_effect
	      REF_LCs[i] = REF_LCs[i]/np.median(REF_LCs[i])
	    master_ref = np.array([np.median(a) for a in zip(*(REF_LCs))])
	    master_ref = master_ref/np.median(master_ref)
	  

	  #bc.bin_calc(folder,N_ref_stars,JD,flux_T,err_T_MAG,flux_target_norm,T_LC,REF_LC,REF_LCs,fwhm,airmass,err_target_w_mean_norm,err_ref_Rs,MAD_factor,ref_number,iteration)
	  
	  c.calc(folder,N_ref_stars,N_bad_refs,airmass,np_airmass,poly2d_AIRMASS,fwhm,JD,T_LC,err_target_w_mean_norm,norm_ref_star_flux_weight_sum_mean,REF_LCs,err_r_norm_sum,err_ref_Rs,master_ref,iteration)
	  iteration += 1
