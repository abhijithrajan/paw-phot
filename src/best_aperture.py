#! /usr/bin/env python 

#=======================================#
#    Written by: Paul Anthony Wilson	#
# 	   paw@astro.ex.ac.uk		#
#=======================================#

import numpy as np
from numpy import loadtxt, array
from os import listdir, getcwd
import sys

from pylab import *

def chi2(obs,err):
  core = ((obs - 1.) / err)**2.
  return core.sum()/len(obs-1)

def best_app(folder,one_or_all):
  if one_or_all == 'on':
    workdir = getcwd()
    aperture_dirs = [x for x in listdir(workdir+'/'+folder) if x.startswith('A') ]
    aperture_dirs.sort()
    ds9 = loadtxt(workdir+'/'+folder+'/data/ds9.reg')
  else:
    workdir = getcwd()
    aperture_dirs = [x for x in listdir(workdir+'/'+folder) if x.startswith('A') ]
    aperture_dirs.sort()
    ds9 = loadtxt(workdir+'/'+folder[:-8]+'data/ds9.reg')
  
  N_ref_stars = len(ds9)-1		# Defining the number of reference star

  ''' Initialising 2D arrays '''
  flux_T = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  flux_T_MAG = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  flux_target_norm = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  err_T = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  err_target_norm = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  continuum = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]

  flux_r = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]		# Creating 2D lists
  mag_r = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  err_r = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  err_r_MAG = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  flux_r_norm = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  err_r_norm = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  stddev = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  combined_ref_stars = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  stddev = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  stddev_norm = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]

  #norm_ref_star_flux_weight_sum = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  flux_no_weight  = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]

  variance = [[] for _ in range(N_ref_stars)]
  variance_norm = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  one_over_variance = [[] for _ in range(N_ref_stars)]

  flux_ref_weight_sum = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  norm_ref_star_flux_sum = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  
  sigma_weighted_mean_ref_stars = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
  err_target_w_mean_norm = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]

  best_combo = []
  best_target = []
  stddev_norm_mean_final = [] #mean of all ref stars
  ''' Initialisation done '''
  
  print "\n"
  ap = []
  chi = []
  chi_sq = []
  for g in range(0,len(aperture_dirs)):
    stddev_norm_mean = []
    x, y, flux, mag, merr, msky, xcenter, ycenter = loadtxt(workdir+'/'+folder+"/"+aperture_dirs[g]+'/xyf.list', unpack=True)

    ### NEW  
    flux_T[g] = flux[::N_ref_stars+1]
    flux_T_MAG[g] = mag[::N_ref_stars+1]
    flux_target_norm[g] = flux_T[g]/flux_T[g].mean()
    err_T[g] = merr[::N_ref_stars+1]*flux_T[g]/1.0857
    err_target_norm[g] = err_T[g]/flux_T[g].mean()
    ### NEW END
    for i in range(N_ref_stars):
          flux_r[g][i] = np.array(flux[i+1::N_ref_stars+1])			# flux_r[ref_star][frame]
          flux_r_norm[g][i] = np.array(flux_r[g][i]/flux_r[g][i].mean())		# flux_r normalised
          err_r[g][i] = np.array(merr[i+1::N_ref_stars+1]*flux_r[g][i]/1.0857)				# flux_r[ref_star][frame]
          err_r_norm[g][i] = np.array(err_r[g][i]/flux_r[g][i].mean()+1e-12)
          stddev[g][i] = flux_r_norm[g][i].std()							# STTDEV of reference stars

    flux_r_norm[g] = array(flux_r_norm[g])
    err_r_norm[g] = array(err_r_norm[g])
    variance_norm[g] = array(err_r_norm[g]*err_r_norm[g])


    weights = []
    error_ref =[]
    #wmean = []
    #wmean = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]
    #total_weight = [[[] for _ in range(N_ref_stars)] for _ in range(len(aperture_dirs))]

    ''' Find the weighted mean of the flux of the reference stars '''
    norm_ref_star_flux_weight_sum=np.zeros(len(flux_r_norm[g][0]))
    sigma_weighted_mean_ref_stars=np.zeros(len(flux_r_norm[g][0]))
    for j in range(len(flux_r_norm[g][0])):
      wmean = 0.
      total_weight = 0.
      for i in range(N_ref_stars):
        wmean += (flux_r_norm[g][i][j] / variance_norm[g][i][j])
        total_weight += (1. / variance_norm[g][i][j])

      norm_ref_star_flux_weight_sum[j] = wmean / total_weight
      sigma_weighted_mean_ref_stars[j] = np.sqrt(1./total_weight)


    ''' FINAL ERROR BARS '''
    err_target_w_mean_norm[g] = np.sqrt(err_target_norm[g]**2+sigma_weighted_mean_ref_stars[g]**2)    

    ''' Computing the weighted mean of the reference flux for each ref star'''
    #norm_ref_star_flux_weight_sum[g] = np.average(flux_r_norm[g], axis=0, weights=weights)
    norm_ref_star_flux_sum[g] = array([sum(a) for a in zip(*(flux_r_norm[g]))])
    norm_ref_star_flux_sum[g] = norm_ref_star_flux_sum[g]/N_ref_stars
    stddev_norm_mean = list(stddev_norm_mean)
    for i in range(N_ref_stars):
      stddev_norm_mean.append((norm_ref_star_flux_sum[g]/norm_ref_star_flux_weight_sum).std())
    stddev_norm_mean = array(stddev_norm_mean)

    stddev_norm_mean_final.append(stddev_norm_mean.mean() )
    
    chi_sq.append(chi2(norm_ref_star_flux_sum[g]/norm_ref_star_flux_weight_sum,sigma_weighted_mean_ref_stars[g]))

    print aperture_dirs[g],"\t-->\t",stddev_norm_mean_final[g],"\t-->\t",chi_sq[g]
    best_combo.append([stddev_norm_mean_final[g],chi_sq[g],aperture_dirs[g]])
    ap.append(aperture_dirs[g][1:5])
    chi.append(chi_sq[g])
    
  best_combo = min(best_combo)
    
  if one_or_all == 'on':
	print "\nBest Ref. aperture: ",best_combo[2]," with a Continuum STDEV = ",best_combo[0]," and a Reduced Chi^2 of ",best_combo[1]
  
	plt.figure(1)
	ax = subplot('211')
	plot(ap,chi,'.k')
	#plot(np.arange(len(chi)),chi,'.k')
	ax2 = subplot('212')
	plot(ap,stddev_norm_mean_final,'.r')
	xlabel('JD [hours]')
	ylabel('Relative Flux')
	#savefig(folder[28:]+'_trend.pdf')
	show()
  else:
    print "\nBest Ref. aperture: ",best_combo[2]," with a Continuum STDEV = ",best_combo[0]," and a Reduced Chi^2 of ",best_combo[1]
    f = open('best_apertures.txt', 'a')
    print >> f, folder[9:-9],best_combo[2]," with a Continuum STDEV = ",best_combo[0]," and a Reduced Chi^2 of ",best_combo[1]
    f.close()
    plt.figure(1)
    ax = subplot('211')
    cla()
    plot(ap,chi,'.k')
    ax2 = subplot('212')
    cla()
    plot(ap,stddev_norm_mean_final,'.r')
    xlabel('JD [hours]')
    ylabel('Relative Flux')
    savefig(folder[9:-9]+'_aperture.pdf')
	#show()