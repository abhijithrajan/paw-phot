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
from os import listdir
import os
import time
import matplotlib
import numpy.ma as ma
from scipy import stats

import math
import random
import unix as u
from scipy.signal import spectral
from itertools import izip_longest

ion()

def get_airmass_JD(folder):
  global JD, imgs
  dir_contents = listdir(folder+'/data')
  imgs = [s for s in dir_contents if s.startswith("ali") and s.endswith(".fits")]
  imgs = sorted(imgs)
  headers = []
  airmass = []
  JD = []
  for i in range(len(imgs)):
    headers.append(pyfits.open(folder+'/data/'+imgs[i]))
    JD.append(headers[i][0].header['MJD-OBS'])
    airmass.append(headers[i][0].header['HIERARCH ESO TEL AIRM START'])
  return airmass, JD

def chi2(obs,err):
  obs = obs/np.median(obs)
  core = ((obs - 1.) / err)**2.
  return core.sum()/len(obs-1)

def robust(obs,err):
  core = abs((obs-median(obs))/err)
  return core.sum()/len(obs-1)

def coloumn_weighted_mean(B,err):
  data = {}
  weights_N = {}
  results = []
  max_len = 0
  
  weights = [[] for _ in range(N_ref_stars)]
  for i in range(len(err)):
    err[i] = np.array(err[i])
    variance = err[i]*err[i]
    weights[i] = 1./variance
    weights[i] = weights[i]/weights[i].sum()

  for alist in B:
	  length = len(alist)
	  max_len = length if (length > max_len) else max_len
	  for i in range(length):
		  data.setdefault(i, []).append(alist[i])

  for alist in weights:
	  length = len(alist)
	  max_len = length if (length > max_len) else max_len
	  for i in range(length):
		  weights_N.setdefault(i, []).append(alist[i])
		  
  for i in range(max_len):
	  vals = np.array(data[i])
	  vals_w = np.array(weights_N[i])
	  results.append(np.average(vals, weights=vals_w))
  
  return np.array(results)

def coloumn_median(B,A):
  maxlen = max(len(x) for x in B)
  if len(A) == maxlen:
    ans = A
  else:
    C = np.array([l+[np.nan]*(maxlen-len(l)) for l in B])
    dat = np.ma.masked_array(C,np.isnan(C))
    ans = np.median(dat,axis=0)
  return ans

def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default
    """
    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m

def ref_test(folder,JD,flux_T,flux_r_norm):
	name = folder[33:39]
	fig=plt.figure(1,figsize=(11.69,8.27))
	f1=fig.add_subplot(1,1,1)
	f1.cla()
	
	# Plot the refs
	for i in range(len(flux_r_norm)):
	  f1.plot(JD,flux_r_norm[i],'.',label='Faint')
	
	plt.title('CONSTANT Aperture')
	plt.xlabel('JD')
	plt.ylabel('Flux')
	plt.ylim(0.8,1.2)
	plt.tight_layout()
	plt.savefig(folder+'/tests/'+name+'_flux_vs_time_CONSTANT_AP.pdf')

def bin_calc(folder,JD,flux_T,err_T_MAG,flux_target_norm,T_LC,REF_LC,REF_LCs,fwhm,airmass,err_target_w_mean_norm,err_ref_Rs,MAD_factor,ref_number,iteration):
	
	global master_entire_ref_LC_detrended
	global flux_entire_detrended
	global flux_binned_detrended
	global flux_entire_Rs_detrended
	global flux_binned_C_Rs_detrended
	
	
	dividers = []	
	out_i = []
	
	airmass_bin = []
	airmass_bin_C = []
	airmass_binned_C = []
	
	fwhm_bin = []
	fwhm_bin_C = []
	fwhm_binned_C = []
	
	airmass_entire = []
	fwhm_entire = []
	
	bad_ref_stars = []	
	good_ref_index = []

	''' Initialising target star arrays '''	
	JD_bin = []
	flux_bin = []
	JD_bin_mean = []
	flux_bin_mean = []
	JD_bin_C = []
	flux_bin_C = []
	
	JD_binned_C= []
	flux_binned_C = []
	flux_binned_C_err = []

	JD_bin_removed = []	
	flux_bin_removed = []

	JD_entire = []	
	flux_entire = []

	JD_binned_removed_entire = []
	flux_binned_removed_entire = []
	
	bin_std = []
	MAD_bin = []
	
	chi2_good_refs = []

	err_target_w_mean_norm_bin = []
	err_target_w_mean_norm_bin_C = []
	err_target_w_mean_norm_binned_C = []

	''' Initialising reference star arrays '''	
	JD_bin_Rs = [[] for _ in range(N_ref_stars)]
	flux_bin_Rs = [[] for _ in range(N_ref_stars)]
	
	fwhm_bin_Rs = [[] for _ in range(N_ref_stars)]
	airmass_bin_Rs = [[] for _ in range(N_ref_stars)]
	
	JD_bin_mean_Rs = [[] for _ in range(N_ref_stars)]
	flux_bin_mean_Rs = [[] for _ in range(N_ref_stars)]

	JD_bin_C_Rs = [[] for _ in range(N_ref_stars)]
	flux_bin_C_Rs = [[] for _ in range(N_ref_stars)]
	
	fwhm_bin_C_Rs = [[] for _ in range(N_ref_stars)]
	airmass_bin_C_Rs = [[] for _ in range(N_ref_stars)]

	JD_bin_removed_Rs = [[] for _ in range(N_ref_stars)]
	flux_bin_removed_Rs = [[] for _ in range(N_ref_stars)]
	
	JD_entire_Rs = [[] for _ in range(N_ref_stars)]
	flux_entire_Rs = [[] for _ in range(N_ref_stars)]
	flux_entire_Rs_err = [[] for _ in range(N_ref_stars)]
	
	fwhm_entire_Rs = [[] for _ in range(N_ref_stars)]
	airmass_entire_Rs = [[] for _ in range(N_ref_stars)]

	JD_binned_C_Rs = [[] for _ in range(N_ref_stars)]	
	flux_binned_C_Rs = [[] for _ in range(N_ref_stars)]
	flux_binned_C_Rs_err = [[] for _ in range(N_ref_stars)]
	
	fwhm_binned_C_Rs =  [[] for _ in range(N_ref_stars)]
	airmass_binned_C_Rs =  [[] for _ in range(N_ref_stars)]
	
	linear_FWHM_entire_Rs =  [[] for _ in range(N_ref_stars)]
	poly2d_AIRMASS_entire_Rs =  [[] for _ in range(N_ref_stars)]
	fitted_trend_entire_Rs =  [[] for _ in range(N_ref_stars)]
	
	flux_entire_Rs_detrended = [[] for _ in range(N_ref_stars)]
	flux_binned_Rs_detrended = [[] for _ in range(N_ref_stars)]
	flux_binned_Rs_detrended_err = [[] for _ in range(N_ref_stars)]
	
	JD_binned_removed_entire_Rs = [[] for _ in range(N_ref_stars)]
	flux_binned_removed_entire_Rs = [[] for _ in range(N_ref_stars)]

	bin_std_Rs = [[] for _ in range(N_ref_stars)]
	MAD_bin_Rs = [[] for _ in range(N_ref_stars)]

	err_ref_bin_Rs = [[] for _ in range(N_ref_stars)]
	err_ref_bin_C_Rs = [[] for _ in range(N_ref_stars)]
	err_ref_binned_C_Rs = [[] for _ in range(N_ref_stars)]
	 
	power_Rs = [[] for _ in range(N_ref_stars)]
	
	correlation = [[] for _ in range(N_ref_stars)]
	
	''' Drawing the bins. Algorithm based on calculating the time difference
		between each exposure. If a jump is detected a new bin line is drawn
		'''
	divide = []
	index_divide = [0]
	for i in range(len(JD)-2):
	  if (JD[i+2]-JD[i+1]) >= (JD[i+1]-JD[i])*4.0:
	    index_divide.append(i+2)
	    divide.append(JD[i+1] + (JD[i+2]-JD[i+1])/2.)

	index_divide.append(i+2)
	divide_at = (JD[-1] - JD[0])/5.
	time_step = np.arange(JD[0],JD[-1],divide_at)
	k = 1
	if len(index_divide) < 3:
	  index_divide = [0]
	  for i in range(len(JD-2)):
	    if (divide_at*k <= JD[i] <= divide_at*(k+1)) and (k < len(time_step)-1):
	      print divide_at*k,JD[i],divide_at*(k+1)
	      index_divide.append(i)
	      divide.append(time_step[k])
	      k += 1	

	''' SPECIAL CASES '''	
	if folder[30:32] == '07':
	  del index_divide[4]
	  index_divide.append(80)
	  index_divide.append(100)
	  index_divide.append(120)
	  index_divide.append(140)
	  index_divide.append(160)
	  index_divide.append(180)
	  index_divide.append(200)
	  index_divide.append(217)
	  divide.append(JD[80+1] + (JD[80+2]-JD[80+1])/2.)
	  divide.append(JD[100+1] + (JD[100+2]-JD[100+1])/2.)
	  divide.append(JD[120+1] + (JD[120+2]-JD[120+1])/2.)
	  divide.append(JD[140+1] + (JD[140+2]-JD[140+1])/2.)
	  divide.append(JD[160+1] + (JD[160+2]-JD[160+1])/2.)
	  divide.append(JD[180+1] + (JD[180+2]-JD[180+1])/2.)
	  divide.append(JD[200+1] + (JD[200+2]-JD[200+1])/2.)
	  divide.append(JD[217+1] + (JD[217+2]-JD[217+1])/2.)

	if folder[30:32] == '12':
	  del index_divide[-1]
	  del divide[-1]

	if folder[30:32] == '17':
	  del index_divide[-1]
	  del divide[-1]

	'''
	if folder[30:32] == '35':
	  del index_divide[12]
	  print JD[300+1] + (JD[300+2]-JD[300+1])/2.
	  print index_divide
	  index_divide.append(300)
	  divide.append(JD[300+1] + (JD[300+2]-JD[300+1])/2.)
	  index_divide.append(322)
	  divide.append(JD[322+1] + (JD[322+2]-JD[322+1])/2.)
	'''
	
	if folder[30:32] == '37':
	  del index_divide[-1]
	  del divide[-1]

	if folder[30:32] == '46':
	  del index_divide[5]
	  del divide[5]
	

	''' The following lines of code treats each data bin individually, figuring out which points are
		outliers. It then calculates a mean for the data within each bin before finally outputting
		the unbinned lightcurve again with outliers removed.'''
	for j in range (len(index_divide)-1):					# Since we start counting at 0
	  outliers_i = []
	  outliers_i_R = []
	  outliers_i_Rs = [[] for _ in range(N_ref_stars)]
	  
	  dividers.append(divide[j-1])						# Used to draw the vertical lines

	  JD_bin.append(np.array(JD[index_divide[j]:index_divide[j+1]]))	# Define which points are within the bin
	  JD_bin_mean.append(JD_bin[j].mean())					# Calculate mean of points within the bin
	  
	  flux_bin.append(T_LC[index_divide[j]:index_divide[j+1]])		# Define which points are within the bin
	  flux_bin_mean.append(flux_bin[j].mean())				# Calculate mean of points within the bin
	  	  
	  bin_std.append(flux_bin[j].std())					# STTDEV of flux points within each bin. Used to calculated errors later.
	  MAD_bin.append(MAD(flux_bin[j], c=0.6745, axis=None))			# The MAD of each bin
	  
	  fwhm_bin.append(fwhm[index_divide[j]:index_divide[j+1]])
	  airmass_bin.append(airmass[index_divide[j]:index_divide[j+1]])
	  
	  err_target_w_mean_norm_bin.append(err_target_w_mean_norm[index_divide[j]:index_divide[j+1]])
	  
	  ''' Same calculation as above but for the reference stars '''
	  for i in range(N_ref_stars):
	    JD_bin_Rs[i].append(np.array(JD[index_divide[j]:index_divide[j+1]]))
	    JD_bin_mean_Rs[i].append(JD_bin_Rs[i][j].mean())
	    flux_bin_Rs[i].append(REF_LCs[i][index_divide[j]:index_divide[j+1]])
	    flux_bin_mean_Rs[i].append(flux_bin_Rs[i][j].mean())
	    
	    fwhm_bin_Rs[i].append(np.array(fwhm[index_divide[j]:index_divide[j+1]]))
	    airmass_bin_Rs[i].append(np.array(airmass[index_divide[j]:index_divide[j+1]]))
	    
	    bin_std_Rs[i].append(flux_bin_Rs[i][j].std())
	    MAD_bin_Rs[i].append(MAD(flux_bin_Rs[i][j], c=0.6745, axis=None))
	    
	    err_ref_bin_Rs[i].append(err_ref_Rs[i][index_divide[j]:index_divide[j+1]])
	    
	  ''' Outliers above the threshold MAD_bin * MAD_factor
	  	  are saved in the outliers list, and later plotted as grey crosses
	  	  '''
	  for k in range(len(flux_bin[j])):
	    
	    # Find outliers in the target light curve. Working with a median as an average would be influenced by outliers.
	    if  (abs(np.median(flux_bin[j])-flux_bin[j][k]) >= MAD_bin[j]*MAD_factor):
	  	  outliers_i.append(k)
	  	
	    # Find outliers in the reference star light curves
	    for i in range(N_ref_stars):
	      if  (abs(np.median(flux_bin_Rs[i][j])-flux_bin_Rs[i][j][k]) >= MAD_bin_Rs[i][j]*MAD_factor):
	  	    outliers_i_Rs[i].append(k)

	  ''' Create lists of the values of the removed outliers'''
	  JD_bin_removed.append(JD_bin[j][outliers_i])
	  flux_bin_removed.append(flux_bin[j][outliers_i])
  
	  for i in range(N_ref_stars):
	    JD_bin_removed_Rs[i].append(JD_bin_Rs[i][j][outliers_i_Rs[i]])
	    flux_bin_removed_Rs[i].append(flux_bin_Rs[i][j][outliers_i_Rs[i]])
	  
	  ''' Creating lists with the outliers removed.
	  	  _C indicates that outliers have been removed '''
	  JD_bin_C.append([x for x in JD_bin[j] if x not in JD_bin[j][outliers_i]])
	  flux_bin_C.append([x for x in flux_bin[j] if x not in flux_bin[j][outliers_i]])

	  fwhm_bin_C.append([x for x in fwhm_bin[j] if x not in fwhm_bin[j][outliers_i]])
	  airmass_bin_C.append([x for x in airmass_bin[j] if x not in airmass_bin[j][outliers_i]])	  
	  
	  for i in range(N_ref_stars):
	    JD_bin_C_Rs[i].append([x for x in JD_bin_Rs[i][j] if x not in JD_bin_Rs[i][j][outliers_i_Rs[i]]])		#JD_bin would be the same as JD_bin_R
	    flux_bin_C_Rs[i].append([x for x in flux_bin_Rs[i][j] if x not in flux_bin_Rs[i][j][outliers_i_Rs[i]]])
	    err_ref_bin_C_Rs[i].append([x for x in  err_ref_bin_Rs[i][j] if x not in  err_ref_bin_Rs[i][j][outliers_i_Rs[i]]])
	    
	    fwhm_bin_C_Rs[i].append([x for x in  fwhm_bin_Rs[i][j] if x not in  fwhm_bin_Rs[i][j][outliers_i_Rs[i]]])
	    airmass_bin_C_Rs[i].append([x for x in  airmass_bin_Rs[i][j] if x not in  airmass_bin_Rs[i][j][outliers_i_Rs[i]]])

	  err_target_w_mean_norm_bin_C.append([x for x in err_target_w_mean_norm_bin[j] if x not in err_target_w_mean_norm_bin[j][outliers_i]])
	  
	  ''' Converting corrected bin data into arrays for further calculations '''
	  JD_bin_C[j] = np.array(JD_bin_C[j])
	  flux_bin_C[j] = np.array(flux_bin_C[j])

	  fwhm_bin_C[j] = np.array(fwhm_bin_C[j])
	  airmass_bin_C[j] = np.array(airmass_bin_C[j])
	  err_target_w_mean_norm_bin_C[j] = np.array(err_target_w_mean_norm_bin_C[j])

	  for i in range(N_ref_stars):
	    JD_bin_C_Rs[i][j] = np.array(JD_bin_C_Rs[i][j])
	    flux_bin_C_Rs[i][j] = np.array(flux_bin_C_Rs[i][j])
	    err_ref_bin_C_Rs[i][j] = np.array(err_ref_bin_C_Rs[i][j])
	    
	    fwhm_bin_C_Rs[i][j] = np.array(fwhm_bin_C_Rs[i][j])
	    airmass_bin_C_Rs[i][j] = np.array(airmass_bin_C_Rs[i][j])
	  	  
	  ''' Calculating the mean of the points within each bin.'''
	  JD_binned_C.append(np.median(JD_bin_C[j]))#JD_bin_C[j].mean())
	  flux_binned_C.append(np.median(flux_bin_C[j]))#flux_bin_C[j].mean())
	  fwhm_binned_C.append(np.median(fwhm_bin_C[j]))#fwhm_bin_C[j].mean())
	  err_target_w_mean_norm_binned_C.append(np.median(err_target_w_mean_norm_bin_C[j]))#/np.sqrt(len(err_target_w_mean_norm_bin_C[j])))#err_target_w_mean_norm_bin_C[j].mean())	# Standard deviation of the mean based upon IRAF errors
	
	  airmass_binned_C.append(np.median(airmass_bin_C[j]))#airmass_bin_C[j].mean())
	  flux_binned_C_err.append(flux_bin_C[j].std())#/np.sqrt(len(flux_bin_C[j])))	# Standard deviation of the mean. Errors from the bins themselves.

	  
	  for i in range(N_ref_stars):
	    JD_binned_C_Rs[i].append(np.median(JD_bin_C_Rs[i][j]))#JD_bin_C_Rs[i][j].mean())
	    flux_binned_C_Rs[i].append(np.median(flux_bin_C_Rs[i][j]))#flux_bin_C_Rs[i][j].mean())
	    flux_binned_C_Rs_err[i].append(flux_bin_C_Rs[i][j].std())#/np.sqrt(len(flux_bin_C_Rs[i][j])))		 # Standard deviation of the mean. Errors from the bins themselves.
	    err_ref_binned_C_Rs[i].append(np.median(err_ref_bin_C_Rs[i][j]))#/np.sqrt(len(err_ref_bin_C_Rs[i][j])))# Standard deviation of the mean. IRAF errors
	    
	    fwhm_binned_C_Rs[i].append(np.median(fwhm_bin_C[j]))#fwhm_bin_C_Rs[i][j].mean())
	    airmass_binned_C_Rs[i].append(np.median(airmass_bin_C_Rs[i][j]))#airmass_bin_C_Rs[i][j].mean())
	  
	  ''' Adding data from the bins together to create continuous unbinned corrected
	  	  light curves designated with the _entire label.
	  	  '''	  	  
	  # For the target
	  for k in range(len(JD_bin_C[j])):
	    JD_entire.append(JD_bin_C[j][k])
	    flux_entire.append(flux_bin_C[j][k])
	    fwhm_entire.append(fwhm_bin_C[j][k])
	    airmass_entire.append(airmass_bin_C[j][k])
	    
	  # Creating a light curve of the removed target points
	    try:
	      if JD_bin_removed[j][k]:
	        JD_binned_removed_entire.append(JD_bin_removed[j][k])
	        flux_binned_removed_entire.append(flux_bin_removed[j][k])
	    except:
	      pass


	  # For the refs
	  for i in range(N_ref_stars):
	    for k in range(len(JD_bin_C_Rs[i][j])):
	      JD_entire_Rs[i].append(JD_bin_C_Rs[i][j][k])
	      flux_entire_Rs[i].append(flux_bin_C_Rs[i][j][k])
	      flux_entire_Rs_err[i].append(err_ref_bin_C_Rs[i][j][k])
	      fwhm_entire_Rs[i].append(fwhm_bin_C_Rs[i][j][k])
	      airmass_entire_Rs[i].append(airmass_bin_C_Rs[i][j][k])

	      # Creating a light curve of the removed reference points
	      try:
	        if JD_bin_removed_Rs[i][j][k]:
	          JD_binned_removed_entire_Rs[i].append(JD_bin_removed_Rs[i][j][k])
	          flux_binned_removed_entire_Rs[i].append(flux_bin_removed_Rs[i][j][k])
	      except:
	        pass

	dividers = np.array(dividers)
	bin_std = np.array(bin_std)

	''' The mean of the reference star light curves and the FWHM/AIRMASS (with points removed).
		This is used for calculating airmass and FWHM trends '''
	
	JD_entire_Rs = np.array(JD_entire_Rs)
	flux_entire_Rs = np.array(flux_entire_Rs)
	flux_entire_Rs_err = np.array(flux_entire_Rs_err)
	fwhm_entire_Rs = np.array(fwhm_entire_Rs)
	airmass_entire_Rs = np.array(airmass_entire_Rs)
	
	for i in range(N_ref_stars):
	  JD_binned_C_Rs[i] = np.array(JD_binned_C_Rs[i])
	  JD_binned_removed_entire_Rs[i] = np.array(JD_binned_removed_entire_Rs[i])
	  flux_binned_removed_entire_Rs[i] = np.array(flux_binned_removed_entire_Rs[i])	  
	  fwhm_entire_Rs[i] = np.array(fwhm_entire_Rs[i])
	  airmass_entire_Rs[i] = np.array(airmass_entire_Rs[i])	
	
	master_entire_ref_JD = coloumn_median(JD_entire_Rs,JD)
	master_entire_ref_LC = coloumn_weighted_mean(flux_entire_Rs,flux_entire_Rs_err)	# Flux_entire_Rs_err is based in IRAF ERRORS
	master_entire_FWHM = np.array([sum(col) / sum(cmp(x,0) for x in col) for col in izip_longest(*fwhm_entire_Rs, fillvalue=0)])
	master_entire_AIRMASS = np.array([sum(col) / sum(cmp(x,0) for x in col) for col in izip_longest(*airmass_entire_Rs, fillvalue=0)])
	
	''' Calculating the linear slopes of Normalised flux vs FWHM and AIRMASS '''
	if (iteration == 0):
	  slope_FWHM, intercept_FWHM, r_value_FWHM, p_value_FWHM, std_err_FWHM = stats.linregress(master_entire_FWHM,master_entire_ref_LC)#(fwhm_entire,flux_entire)#
	  slope_AIRMASS, intercept_AIRMASS, r_value_AIRMASS, p_value_AIRMASS, std_err_AIRMASS = stats.linregress(master_entire_AIRMASS,master_entire_ref_LC)
	  A,B,C = np.polyfit(master_entire_AIRMASS,master_entire_ref_LC, 2) # Fitting a 2D polynomial to the airmass
	else:
	  slope_FWHM, intercept_FWHM, r_value_FWHM, p_value_FWHM, std_err_FWHM = stats.linregress(master_entire_FWHM,master_entire_ref_LC_detrended)
	  slope_AIRMASS, intercept_AIRMASS, r_value_AIRMASS, p_value_AIRMASS, std_err_AIRMASS = stats.linregress(master_entire_AIRMASS,master_entire_ref_LC_detrended)
	  A,B,C = np.polyfit(master_entire_AIRMASS,master_entire_ref_LC_detrended, 2) # Fitting a 2D polynomial to the airmass
	
	''' Converting from list to arrays for further calculations '''
	JD_entire = np.array(JD_entire)
	flux_entire = np.array(flux_entire)

	JD_binned_removed_entire = np.array(JD_binned_removed_entire)
	flux_binned_removed_entire = np.array(flux_binned_removed_entire)
	
	fwhm_entire = np.array(fwhm_entire)
	airmass_entire = np.array(airmass_entire)

	fwhm_binned_C = np.array(fwhm_binned_C)
	airmass_binned_C = np.array(airmass_binned_C)
	err_target_w_mean_norm_binned_C = np.array(err_target_w_mean_norm_binned_C)
	
	for i in range(N_ref_stars):
	  flux_entire_Rs[i] = np.array(flux_entire_Rs[i])
	  err_ref_binned_C_Rs[i] = np.array(err_ref_binned_C_Rs[i])
	  linear_FWHM_entire_Rs[i] = slope_FWHM*fwhm_entire_Rs[i]+intercept_FWHM
	  poly2d_AIRMASS_entire_Rs[i] = A*airmass_entire_Rs[i]**2+B*airmass_entire_Rs[i]+C
	
	''' Defining linear FWHM trend '''
	linear_FWHM = slope_FWHM*master_entire_FWHM+intercept_FWHM
	linear_FWHM_entire = slope_FWHM*fwhm_entire+intercept_FWHM
	linear_binned_FWHM = slope_FWHM*fwhm_binned_C+intercept_FWHM
	
	''' Defining linear AIRMASS trend '''
	linear_AIRMASS = slope_AIRMASS*master_entire_AIRMASS+intercept_AIRMASS
	linear_AIRMASS_entire = slope_AIRMASS*airmass_entire+intercept_AIRMASS
	
	''' Calculating 2D polynomial AIRMASS trend '''
	poly2d_AIRMASS = A*master_entire_AIRMASS**2+B*master_entire_AIRMASS+C
	poly2d_AIRMASS_entire = A*airmass_entire**2+B*airmass_entire+C
	poly2d_binned_AIRMASS = A*airmass_binned_C**2+B*airmass_binned_C+C
	
		
	''' Calculating the difference in Y-axis to determine which
		linear trend is dominating. '''
	y_diff_AIRMASS = linear_AIRMASS.max()-linear_AIRMASS.min()
	y_diff_FWHM = linear_FWHM.max()-linear_FWHM.min()

	if abs(y_diff_FWHM) >= abs(y_diff_AIRMASS):
	  print "\nFWHM trends dominate ","(",round(y_diff_FWHM/y_diff_AIRMASS,1),")"  
	  dominating_trend = 'FWHM'
	  fitted_trend = linear_FWHM
	  fitted_trend_entire = linear_FWHM_entire
	  fitted_trend_binned = linear_binned_FWHM
	  fitted_trend_secondary = poly2d_AIRMASS
	  for i in range(N_ref_stars):
	    fitted_trend_entire_Rs[i] = linear_FWHM_entire_Rs[i]
	else:
	  print "\nAirmass trends dominate","(",round(y_diff_AIRMASS/y_diff_FWHM,1),")"
	  dominating_trend = 'Airmass'
	  fitted_trend = poly2d_AIRMASS
	  fitted_trend_entire = poly2d_AIRMASS_entire
	  fitted_trend_binned = poly2d_binned_AIRMASS
	  fitted_trend_secondary = linear_FWHM
	  for i in range(N_ref_stars):
	    fitted_trend_entire_Rs[i] = poly2d_AIRMASS_entire_Rs[i]
		
	''' The final binned target light curve detrended '''
	for i in range(N_ref_stars):
	    flux_entire_Rs_detrended[i] = flux_entire_Rs[i]/fitted_trend_entire_Rs[i]
	    flux_binned_Rs_detrended[i] = flux_binned_C_Rs[i]/fitted_trend_binned
	if (iteration == 0):
	  master_entire_ref_LC_detrended = master_entire_ref_LC/fitted_trend
	  flux_entire_detrended = flux_entire/fitted_trend_entire
	  flux_binned_detrended = flux_binned_C/fitted_trend_binned
	else:
	  master_entire_ref_LC_detrended = master_entire_ref_LC_detrended/fitted_trend
	  flux_entire_detrended = flux_entire_detrended/fitted_trend_entire
	  flux_binned_detrended = flux_binned_detrended/fitted_trend_binned
	  
	  for i in range(N_ref_stars):
	    flux_entire_Rs_detrended[i] = flux_entire_Rs_detrended[i]/fitted_trend_entire_Rs[i]
	    flux_entire_Rs_detrended[i] = flux_entire_Rs_detrended[i]/np.median(flux_entire_Rs_detrended[i])
	    flux_binned_Rs_detrended[i] = flux_binned_Rs_detrended[i]/fitted_trend_binned
	    flux_binned_Rs_detrended[i] = flux_binned_Rs_detrended[i]/np.median(flux_binned_Rs_detrended[i])


	''' Normalising LC after detrending '''
	flux_binned_removed_entire = flux_binned_removed_entire/np.median(flux_binned_removed_entire)
	flux_entire_detrended = flux_entire_detrended/np.median(flux_entire_detrended)
	flux_binned_detrended = flux_binned_detrended/np.median(flux_binned_detrended)
	
	flux_binned_Rs_detrended_err = np.array([np.std(a) for a in zip(*(flux_binned_Rs_detrended))])
	
	std_refs = []
	MAD_refs = []
	sort_refs = []
	for i in range(N_ref_stars):
	  std_refs.append(flux_entire_Rs_detrended[i].std())
	  MAD_refs.append(MAD(flux_entire_Rs_detrended[i], c=0.6745, axis=None))
	  sort_refs.append([std_refs[i],i])
	std_refs = sorted(std_refs)
	std_refs_sorted = sorted(sort_refs)
	
	''' CHECK THIS '''
	if (iteration == 1):
	  bad_ref_stars = np.loadtxt('temp.txt',usecols=(0,),ndmin=1)
	#'''
	  for i in range(len(bad_ref_stars)):
	    print "Removing contribution of reference star with index: ",int(bad_ref_stars[i])
	    del std_refs_sorted[-1]
	#'''
	''' CHECK THIS '''
	
	worst_ref = std_refs_sorted[-1][1]
	
	std_refs = np.array(std_refs)
	MAD_refs = np.array(MAD_refs)
	ave_std_refs = std_refs.mean()
	std_std_refs = std_refs.std()
	ave_MAD_refs = MAD_refs.mean()
	std_MAD_refs = MAD_refs.std()
	
	flux_binned_C_err = np.array(flux_binned_C_err)	

	''' Choosing the largest error between the std of the bin derived error and the IRAF error '''
	for i in range(len(flux_binned_C_err)):
	  flux_binned_C_err[i] = max([flux_binned_C_err[i],err_target_w_mean_norm_binned_C[i]])
	  flux_binned_C_Rs_err[ref_number][i] = max(flux_binned_C_Rs_err[ref_number][i],err_ref_binned_C_Rs[ref_number][i])
	  flux_binned_C_Rs_err[worst_ref][i] = max(flux_binned_C_Rs_err[worst_ref][i],err_ref_binned_C_Rs[worst_ref][i])
	
	chi2_target_binned = chi2(flux_binned_detrended,flux_binned_C_err)
	chi2_target_IRAF_binned = chi2(flux_binned_detrended,err_target_w_mean_norm_binned_C)
	chi2_ref_binned = chi2(flux_binned_Rs_detrended[ref_number],flux_binned_C_Rs_err[ref_number])
	chi2_worst_binned = chi2(flux_binned_Rs_detrended[worst_ref],flux_binned_C_Rs_err[worst_ref])
	DOF = len(flux_binned_detrended)-1
	
	robust_target = robust(flux_binned_detrended,flux_binned_C_err)
	robust_sample = robust(flux_binned_Rs_detrended[ref_number],flux_binned_C_Rs_err[ref_number])
		  
	print "\nTarget:"
	print "================================"
	print "Average Err (Binned):\t",round(np.array(flux_binned_C_err).mean(),5)
	print "Chi^2_red (Binned):\t",round(chi2_target_binned,2)
	print "Robust:\t\t\t",round(robust_target,2)
	print "DOF:\t\t\t",DOF
	print "================================"
	
	print "\nSample Ref with index ",ref_number,":"
	print "================================"
	print "Average Err (Binned):\t",round(np.array(flux_binned_C_Rs_err[ref_number]).mean(),5)
	print "Chi^2_red (Binned):\t",round(chi2_ref_binned,2)
	print "Robust:\t\t\t",round(robust_target,2)
	print "DOF:\t\t\t",len(flux_binned_Rs_detrended[ref_number])-1
	print "Worst Ref star:\t\t",worst_ref
	print "================================"
	
	print "\nChance this is a variable [%] ",round(100-(1-stats.chi2.cdf(chi2(flux_binned_detrended,flux_binned_C_err)*DOF, DOF))*100.,3),"\n"

	''' This next part calculates the correlations
	between refs stars and the target.'''
	for i in range(N_ref_stars):
	  correlation[i] = stats.pearsonr(T_LC,REF_LCs[i])[0]
	  if abs(correlation[i]) >= 0.5:
		print u.red("Ref: "),i,"\t",correlation[i]
	  else:
		print "Ref: ",i,"\t",correlation[i]
	correlation_max = max(correlation)
	
	print "\nBinned Refs:"
	print "================================================================================================"
	f = open('temp.txt', 'w+')
	#print >> f, 0
	
	std_refs_sorted.sort(key=lambda x: x[1])
	MAD_factor_2 = 2.0
	if len(std_refs) <= 4:	# Min of four stars are required otherwise no ref stars are removed
	  cut_off = 1.e8
	  print "No cut off, due to few stars\n"
	else:
	  cut_off = ave_MAD_refs+MAD_factor_2*std_MAD_refs
	  print "Cut off: ",cut_off," = ",ave_MAD_refs," + ",MAD_factor_2," *",std_MAD_refs
	for i in range(len(std_refs_sorted)):
	  good_ref_index.append(std_refs_sorted[i][1])
	  if std_refs_sorted[i][0] > cut_off:
	    print >> f, i
	    print u.red("Ref: "),good_ref_index[i],u.red("\tSTDDEV: "),std_refs_sorted[i][0],"\tMAD: ",MAD_refs[i],"\t",u.red("\tChi2_red: "),round(chi2(flux_binned_Rs_detrended[good_ref_index[i]],flux_binned_C_Rs_err[good_ref_index[i]]),2)
	  else:
	    print "Ref: ",good_ref_index[i],"\tSTDDEV: ",std_refs_sorted[i][0],"\tMAD: ",MAD_refs[i],"\t","\tChi2_red: ",round(chi2(flux_binned_Rs_detrended[good_ref_index[i]],flux_binned_C_Rs_err[good_ref_index[i]]),2)
	    chi2_good_refs.append(chi2(flux_binned_Rs_detrended[good_ref_index[i]],flux_binned_C_Rs_err[good_ref_index[i]]))
	f.close()
	chi2_red_worst_good_ref = sorted(chi2_good_refs)[-1]
	
	print "================================================================================================\n"
	if (iteration == 1):
	  print "Removed ",N_ref_stars-len(std_refs_sorted)," of ",N_ref_stars," reference stars (",round(((float(N_ref_stars)-len(std_refs_sorted))/N_ref_stars)*100.,2),"%)\n"

	''' This next part is for calculating periods. Commented out to
	save computation time '''
	'''
	frequencies = np.linspace(2*np.pi/4.0, 2*np.pi/0.8, 1e3)	# Frequencies to be probed 0.8 - 4.0 hours
	power_target_detrended = spectral.lombscargle(JD_entire*24, flux_entire_detrended, frequencies)
	power_ref_detrended = spectral.lombscargle(master_entire_ref_JD*24, master_entire_ref_LC_detrended, frequencies)
	periods = 1./(frequencies / 2.0 / np.pi)
	print "Probing periods in the range: ",periods.min()," to ",periods.max()," hours"

	period_target = 1./(frequencies[np.argmax(power_target_detrended)] / 2.0 / np.pi)
	period_ref = 1./(frequencies[np.argmax(power_ref_detrended)] / 2.0 / np.pi)
	print "Period Target (Detrended)= " + str(period_target),"\n"
	print "Period Comb. Ref. (Detrended) = " + str(1./(frequencies[np.argmax(power_ref_detrended)] / 2.0 / np.pi)),"\n\n"

	for i in range(N_ref_stars):
	  power_Rs[i] = spectral.lombscargle(np.array(JD_entire_Rs[i])*24, flux_entire_Rs_detrended[i] , frequencies)
	  print "Period (Detrended) Ref ",i," = " + str(1./(frequencies[np.argmax(power_Rs[i])] / 2.0 / np.pi))
	'''
	
	phot_err = np.median(bin_std)
	
	variable = u.red("\tScatter too large") if (phot_err >= flux_entire_detrended.std()) else u.green("\tPossibly variable")
	good_master = u.red("\tPoor Master Ref.") if (master_entire_ref_LC_detrended.mean() <= flux_entire_detrended.std()) else u.green("\tGood Master Ref")
	
	variable_print = "Scatter_too_large" if (phot_err >= flux_entire_detrended.std()) else "Possibly_variable"
	good_master_print = "Poor_Master_Ref." if (master_entire_ref_LC_detrended.mean() <= flux_entire_detrended.std()) else "Good_Master_Ref"
	
	p2p = (flux_binned_detrended.max()-flux_binned_detrended.min())*100.
	
	#master_ref_LC2 = np.array([mean(a) for a in zip(*(flux_binned_Rs_detrended))])		# Consider doing a weighted mean? Done.
	master_ref_LC = coloumn_weighted_mean(flux_binned_Rs_detrended,flux_binned_C_Rs_err)
	
	CDF = stats.chi2.cdf((chi2_target_binned-chi2_ref_binned), DOF)	# Actually 1 - CDF (Cumulative Distribution Function)
	#print "HEREERE",master_ref_LC.std()
	#sys.exit()
	
	####
	print "Summary:\n"
	print "\nPhotometric quality:\t\t\t",round(phot_err,4)
	print "Stddev of unbinned Target LC:\t\t",round(flux_entire_detrended.std(),4)
	print "Stddev of LC - Photmetric Quality:\t",round(flux_entire_detrended.std()-phot_err,4),variable
	print "Photometric quality Master Ref:\t\t",round(master_entire_ref_LC_detrended.std(),4),'\t',good_master
	print "Correlation:\t\t\t\t",round(correlation_max,2),"\n"
	print "Calculated using the binned LCs:"
	print "--------------------------------"
	print "Peak to peak Amplitude:\t\t",round(p2p,2),"%"
	print "Photometric quality:\t\t",round(flux_binned_C_err.mean(),4)
	print "Photometric quality IRAF:\t\t",round(err_target_w_mean_norm_binned_C.mean(),4)
	print "Standard dev. of target LC:\t",round(flux_binned_detrended.std(),4)
	print "Target Chi2_red:\t\t",chi2_target_binned
	print "Target Chi2_red w. IRAF errors:",chi2_target_IRAF_binned
	print "Master ref Chi2_red:\t\t",chi2(master_ref_LC,flux_binned_Rs_detrended_err)
	print "Master ref Chi2_red w. IRAF errors:",chi2(master_ref_LC,err_ref_binned_C_Rs)
	print "Ref Chi2_red:\t\t\t",chi2_ref_binned
	print "Worst Chi2_red:\t\t",chi2_worst_binned
	print "Worst Chi2_red value amongst good refs",chi2_red_worst_good_ref
	print "\nChance of variability (compared to sample ref) [%] ",round(CDF*100.,3),"\n"
	print "Working with: ",folder[27:]
	print 
	####
	
	scaling = np.median(bin_std)/np.median(flux_binned_C_err)
	if (scaling >= 1.):
	  print "Increasing Error Bars"
	else:
	  print "Decreasing Error Bars"
	  
	print "\nscaling:\t",scaling
	print "\nRescaled Chi2_red:\t",chi2(flux_binned_detrended,scaling*flux_binned_C_err)
	
	#for i in range(len(fwhm_entire)):
	#  print flux_entire_detrended[i]#JD_entire[i]*24#fwhm_entire[i]#JD_entire*24,flux_entire_detrended
	
	if (iteration == 1):
	  f = open('data.txt', 'a+')
	  print >> f,folder[27:],flux_T.mean(),np.median(flux_T),err_T_MAG.mean(),np.median(err_T_MAG),1./err_T_MAG.mean(),1./np.median(err_T_MAG),flux_entire_detrended.std(),phot_err,flux_entire_detrended.std()-phot_err,master_entire_ref_LC_detrended.std(),np.median(flux_binned_C_err),correlation_max,chi2_target_binned,chi2(master_ref_LC,flux_binned_Rs_detrended_err),chi2(flux_binned_detrended,scaling*flux_binned_C_err),chi2(flux_binned_detrended/master_ref_LC,np.sqrt(flux_binned_C_err**2+flux_binned_Rs_detrended_err**2)),robust_target,chi2_ref_binned,chi2_worst_binned,chi2_red_worst_good_ref,(chi2_target_binned-chi2_red_worst_good_ref),CDF*100.,DOF,len(JD),len(JD_entire),((len(JD)-len(JD_entire))/float(len(JD)))*100.,dominating_trend,'no_period_calculated',good_master_print,variable_print,p2p
	  f.close()
	  
	  
	  f = open('lc_data/'+folder[27:]+'.txt', 'w+')
	  for i in range(len(JD_binned_C)):
	    print >> f, JD_binned_C[i]*24.,flux_binned_detrended[i],flux_binned_C_err[i],JD_binned_C_Rs[ref_number][i]*24.,flux_binned_Rs_detrended[ref_number][i],flux_binned_C_Rs_err[ref_number][i],master_ref_LC[i],flux_binned_Rs_detrended_err[i]
	  f.close()
	  
	  f = open('lc_data/'+folder[27:]+'_diff.txt', 'w+')
	  for i in range(len(JD_binned_C)):
	    print >> f, flux_binned_detrended[i]/master_ref_LC[i],np.sqrt(flux_binned_C_err[i]**2+flux_binned_Rs_detrended_err[i]**2)#,scaling*flux_binned_C_err[i],JD_binned_C_Rs[ref_number][i]*24.,flux_binned_Rs_detrended[ref_number][i],flux_binned_C_Rs_err[ref_number][i],master_ref_LC[i],flux_binned_Rs_detrended_err[i]
	  f.close()

	flux_binned_err_Rs = np.array(flux_binned_C_Rs_err)

	

	JD_binned_C = np.array([value for value in JD_binned_C if not math.isnan(value)])

	fig=plt.figure(1,figsize=(11.69,8.27))
	#fig.add_subplot(3,1,1)
	top_panel=fig.add_subplot(3,1,1)
	top_panel.cla()
	for i in range(len(dividers)):
	  top_panel.axvline(dividers[i],linewidth=2, ls='--',color='r',alpha=0.1)							# Divider lines
	if (iteration == 0):
	  top_panel.plot(JD_binned_removed_entire,flux_binned_removed_entire,'x',color="black",alpha=0.5)	# Removed points
	  top_panel.plot(JD,T_LC,'.',color="black",alpha=0.2)									# Non-detrended light curve
	#top_panel.plot(JD_entire,flux_entire,'.',color="black",alpha=0.2)									# Non-detrended light curve
	top_panel.plot(JD_entire,flux_entire_detrended,'.',color="red",alpha=1.0)							# Detrended light curve
	top_panel.errorbar(JD_binned_C,flux_binned_detrended,yerr=flux_binned_C_err,fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color="black", ecolor="black")
	#ylim(flux_entire.min(),flux_entire.max())
	ylim(flux_binned_removed_entire.min(),flux_binned_removed_entire.max())
	xlabel('Time [days]', fontsize=14)
	ylabel('Normalised Flux', fontsize=14)
	
	middle_panel=fig.add_subplot(3,1,2, sharex=top_panel, sharey=top_panel)
	middle_panel.cla()
	for i in range(len(std_refs_sorted)):
	  #middle_panel.plot(JD_entire_Rs[i],flux_entire_Rs[i],'.',color="black",alpha=0.05)
	  middle_panel.plot(JD_binned_C_Rs[good_ref_index[i]],flux_binned_Rs_detrended[good_ref_index[i]],'.',alpha=0.2)
	middle_panel.errorbar(JD_binned_C_Rs[ref_number],master_ref_LC,yerr=flux_binned_Rs_detrended_err,fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color="red", ecolor="black")
	#middle_panel.plot(master_entire_ref_JD,master_entire_ref_LC,'*',alpha=0.8)
	#middle_panel.plot(JD,flux_target_norm,'-r')	# Not working for some reason
	#for i in range(N_ref_stars):
	#  middle_panel.plot(np.array(JD_entire_Rs[i])*24, flux_entire_Rs[i],'.')
	ylim(flux_entire.min(),flux_entire.max())
	xlabel('Time [days]', fontsize=14)
	ylabel('Normalised Flux', fontsize=14)

	
	bottom_panel=fig.add_subplot(3,3,7)
	bottom_panel.cla()
	#plot(JD*10.,airmass,'--g')
	if abs(y_diff_FWHM) >= abs(y_diff_AIRMASS):
	  bottom_panel.plot(master_entire_FWHM,master_entire_ref_LC,'.k',alpha=0.2)
	  bottom_panel.plot(master_entire_FWHM,fitted_trend,'-k',alpha=0.2)
	  bottom_panel.plot(master_entire_FWHM,master_entire_ref_LC_detrended,'.r')
	  #bottom_panel.plot(fwhm_entire,fitted_trend,'-r')
	  ylabel('Normalised Flux', fontsize=14)
	  xlabel('FWHM', fontsize=14)	
	else:
	  bottom_panel.plot(master_entire_AIRMASS,master_entire_ref_LC,'.k',alpha=0.3)
	  bottom_panel.plot(master_entire_AIRMASS,fitted_trend,'-k',alpha=0.3)
	  bottom_panel.plot(master_entire_AIRMASS,master_entire_ref_LC_detrended,'.r')
	  #bottom_panel.plot(airmass_entire,linear_AIRMASS,'-b')
	  ylabel('Normalised Flux', fontsize=14)
	  xlabel('Airmass', fontsize=14)	

	bottom_panel=fig.add_subplot(3,3,8,sharey=bottom_panel)
	bottom_panel.cla()
	#plot(JD*10.,airmass,'--g')
	if abs(y_diff_FWHM) <= abs(y_diff_AIRMASS):
	  if (iteration == 0):
	    bottom_panel.plot(master_entire_FWHM,master_entire_ref_LC,'.k',alpha=0.3)
	    bottom_panel.plot(master_entire_FWHM,fitted_trend_secondary,'-k',alpha=0.3)
	    bottom_panel.plot(master_entire_FWHM,master_entire_ref_LC_detrended,'.r')
	  else:
	    bottom_panel.plot(master_entire_FWHM,master_entire_ref_LC_detrended,'.r')
	  xlabel('FWHM', fontsize=14)	
	else:
	  if (iteration == 0):
	    bottom_panel.plot(master_entire_AIRMASS,master_entire_ref_LC,'.k',alpha=0.3)
	    bottom_panel.plot(master_entire_AIRMASS,fitted_trend_secondary,'-k',alpha=0.3)
	    bottom_panel.plot(master_entire_AIRMASS,master_entire_ref_LC_detrended,'.r')
	  else:
	    bottom_panel.plot(master_entire_AIRMASS,master_entire_ref_LC_detrended,'.r')
	  xlabel('Airmass', fontsize=14)	 	
	bottom_panel=fig.add_subplot(3,3,9)
	bottom_panel.cla()
	bottom_panel.plot(JD*24.,fwhm/np.median(fwhm),'ob',alpha=0.5)
	bottom_panel.plot(JD*24.,airmass,'-g',alpha=0.5)
	xlabel('JD (Hours)', fontsize=14)
	ylabel('FWHM (Norm.)/Airmass', fontsize=14)
	
	tight_layout()

	savefig('plots/B/'+folder[27:]+'_1_'+str(iteration)+'_plot.pdf')
	
	fig2=plt.figure(2,figsize=(11.69,8.27))
	top_panel2=fig2.add_subplot(2,1,1)
	top_panel2.cla()
	top_panel2.errorbar(JD_binned_C*24.,flux_binned_detrended,yerr=flux_binned_C_err,fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color="black", ecolor="black") # Binned light curve
	top_panel2.errorbar(JD_binned_C*24.,flux_binned_detrended,yerr=err_target_w_mean_norm_binned_C,fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color="black", ecolor="red") # Binned light curve with IRAF errors
	top_panel2.minorticks_on()
	xlabel('Time [Hours]', fontsize=14)
	ylabel('Normalised Target. Flux', fontsize=14)
	
	middle_panel2=fig2.add_subplot(2,1,2, sharex=top_panel2, sharey=top_panel2)
	middle_panel2.cla()
	if (iteration == 0):
	  middle_panel2.errorbar(JD_binned_C_Rs[worst_ref]*24.,flux_binned_Rs_detrended[worst_ref],yerr=flux_binned_C_Rs_err[worst_ref],fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color="#808080", ecolor="black")
	middle_panel2.errorbar(JD_binned_C_Rs[ref_number]*24.,flux_binned_Rs_detrended[ref_number],yerr=flux_binned_C_Rs_err[ref_number],fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color="red", ecolor="black")
	middle_panel2.minorticks_on() 
	xlabel('Time [Hours]', fontsize=14)
	ylabel('Normalised Ref. Flux', fontsize=14)
	tight_layout()
	savefig('plots/B/'+folder[27:]+'_2_'+str(iteration)+'_plot.pdf')
	
	fig3=plt.figure(3,figsize=(8.27,11.69))
	top_panel3 = fig3.add_subplot(2,2,1)
	top_panel3.cla()
	rainbow_colors = iter(cm.rainbow(np.linspace(0, 1, N_ref_stars)))
	top_panel3.plot(master_entire_ref_JD*24.,master_entire_ref_LC_detrended,'or',alpha=0.5)
	top_panel3.plot(JD_entire*24,flux_entire_detrended,'ok')
	#print len(JD),len(master_entire_ref_JD)
	xlabel('Time [Hours]', fontsize=14)
	ylim(0.95,1.05)

	for i in range(N_ref_stars):
	  middle_panel3 = fig3.add_subplot(N_ref_stars,2,2*(i+1),sharex=top_panel3, sharey=top_panel3)
	  middle_panel3.cla()  
	  JD_entire_Rs[i] = np.array(JD_entire_Rs[i])
	  if (iteration == 0):
	    plot(JD_entire_Rs[i]*24, flux_entire_Rs[i],'.',color=next(rainbow_colors),markeredgecolor='black')
	  else:
	    plot(JD_entire_Rs[i]*24, flux_entire_Rs_detrended[i],'.',color=next(rainbow_colors),markeredgecolor='black')
	    if i in bad_ref_stars:
	      plot(JD_entire_Rs[i]*24, flux_entire_Rs_detrended[i],'.',color="grey",markeredgecolor='black')
	      
	  text(JD_entire_Rs[-1][-1]*24+0.2,np.median(flux_entire_Rs[i]),str(i))
	  minorticks_on() 
   
	xlabel('Time [Hours]', fontsize=14)
	ylim(0.95,1.05)

	bottom_panel3 = fig3.add_subplot(2,2,3)
	bottom_panel3.cla()
	bottom_panel3.plot(JD*24.,fwhm/np.median(fwhm),'ob',alpha=0.5)
	bottom_panel3.plot(JD*24.,airmass,'-g',alpha=0.5)
	xlabel('JD (Hours)', fontsize=14)
	ylabel('FWHM (Norm.)/Airmass', fontsize=14)
	'''
	rainbow_colors = iter(cm.rainbow(np.linspace(0, 1, N_ref_stars)))
	for i in range(N_ref_stars):
	  bottom_panel3.plot(periods,power_Rs[i],'-',color=next(rainbow_colors))
	bottom_panel3.plot(periods,power_ref_detrended,'-r',linewidth=3.0)
	bottom_panel3.plot(periods,power_target_detrended,'-k',linewidth=3.0)
	'''
	#xlabel('Period [Hours]', fontsize=14)
#	'ylabel('Power', fontsize=14)
	#tight_layout()
	savefig('plots/B/'+folder[27:]+'_3_'+str(iteration)+'_plot.pdf')
	draw()
	clf()
	
def do_plot(folder):

	global N_ref_stars
	
	''' Defining the number of reference stars '''
	ds9 = loadtxt(folder+'/data/ds9_'+folder[33:39]+'.reg')
	N_ref_stars = len(ds9)-1

	norm_ref_star_flux_weight_sum_mean_Rs = []
	REF_LCs = []
	err_ref_Rs = []
	mag_refs = []
	bad_mag_refs = []

	MAD_factor = 2.00	#cut off in the bins

	''' Read Julian Date and airmass info from .fits headers '''
	airmass, JD = get_airmass_JD(folder)
	airmass = np.array(airmass)
	JD = np.array(JD)
	fwhm = loadtxt(folder+'/data/fwhm.txt',usecols=(0,))
	turbulence_fwhm = np.ones(len(fwhm)) * np.array([random.random() for _ in xrange(len(fwhm))])*1e-12
	turbulence_airmass = np.ones(len(airmass)) * np.array([random.random() for _ in xrange(len(airmass))])*1e-12
	fwhm = fwhm+turbulence_fwhm		# A tiny pertubation is added as some of the FWHM values measured by IRAF are identical causing more FWHM elements to be revomved than what is desired.
	airmass = airmass+turbulence_airmass

	''' Calculate the seeing '''
	pixel_scale = 0.288#par.pixel_scale[0]
	seeing = fwhm*pixel_scale

	''' Establishing directories '''
	dir_contents = listdir(folder)
	folders = [s for s in dir_contents if s.startswith("B")]
	data_folder = [s for s in dir_contents if s.startswith("data")]
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

	while (iteration < 2):
	  print "\n\nITERATION: ",iteration
	  
	  if len(flux_T) <= 30.:
	    print "Too few points: ",len(flux_T),"\n"
	    break

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
	  
	  if (iteration == 1):
	    bad_ref_stars = np.loadtxt('temp.txt',usecols=(0,),ndmin=1)
	    for i in range(len(bad_ref_stars)):
	      bad_ref = int(bad_ref_stars[i])
	      print u.red("Removing contribution from ref star with index: "),bad_ref
	      err_r_norm[bad_ref]=10000.*err_r_norm[bad_ref]
	      variance_norm[bad_ref]=10000.*variance_norm[bad_ref]

	  ''' Find the mean of the raw flux of the reference stars '''
	  norm_ref_star_flux_sum = np.array([sum(a) for a in zip(*(flux_r_norm))])
	  norm_ref_star_flux_sum = norm_ref_star_flux_sum/N_ref_stars

	  ''' Finding a reference star of similar brightness for comparison sake '''
	  if (iteration == 0):
		for i in range(len(mag_r)):			# Calculating the mean magnitude of each reference star
		  mag_refs.append(np.median(mag_r[i]))#mag_r[i].mean())	# throughout the observing sequence.

		for i in range(len(mag_r)):			# Finding the REF star closest in brightness.
		  if mag_refs[i] == min(mag_refs, key=lambda x:abs(x-flux_T_MAG.mean())):
			ref_number = i					# ref_number indicates the choosen ref star
	  else:
	    bad_ref_stars = np.loadtxt('temp.txt',usecols=(0,),ndmin=1)
	    for i in range(len(bad_ref_stars)):
	      bad_ref = int(bad_ref_stars[i])
	      bad_mag_refs.append(mag_refs[bad_ref])
	    mag_refs = [x for x in mag_refs if x not in bad_mag_refs]

	    for i in range(len(mag_refs)):			# Finding the REF star closest in brightness.
		  if mag_refs[i] == min(mag_refs, key=lambda x:abs(x-np.median(flux_T_MAG))):
		    ref_number = i					# ref_number indicates the choosen ref star

	  print "\nTarget and Comparison star mag difference:\t",round(abs(np.median(flux_T_MAG) - mag_refs[ref_number]),2),'\t(Target Brighter)' if ((np.median(flux_T_MAG) - mag_refs[ref_number]) >= 0) else '\t(Ref Brighter)\n'

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
		weights_Rs_sum = np.array([np.median(a) for a in zip(*(weights_Rs))])
		weighted_mean_Rs = [x for x in weighted_mean_Rs if x not in weighted_mean_Rs[i]] # Removing the contribution by the comparison star
		weighted_mean_Rs_sum = np.array([np.median(a) for a in zip(*(weighted_mean_Rs))])
		norm_ref_star_flux_weight_sum_mean_Rs.append(weighted_mean_Rs_sum / weights_Rs_sum)
		REF_LCs.append(flux_r_norm[i]/norm_ref_star_flux_weight_sum_mean_Rs[i])
		err_ref_Rs.append(np.sqrt(err_r_norm[i]**2 + np.median(sigma_weighted_mean_ref_stars)**2))	# IRAF errors
		
	  
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
	
	  ''' The non detrended but weighted sample reference star light curve '''
	  REF_LC = REF_LCs[ref_number]
	
	  #bin_calc(folder,JD,flux_T,err_T_MAG,flux_target_norm,T_LC,REF_LC,REF_LCs,fwhm,airmass,err_target_w_mean_norm,err_ref_Rs,MAD_factor,ref_number,iteration)
	  
	  #ref_test(folder,JD,flux_T,flux_r_norm)
	  
	  iteration += 1
