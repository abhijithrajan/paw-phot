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
import time
import matplotlib
import numpy.ma as ma
from scipy import stats

import unix as u

ion()

from matplotlib.dates import date2num
from datetime import datetime

def get_airmass_JD(folder):
  global JD, imgs
  dir_contents = listdir(folder+'/data')
  imgs = [s for s in dir_contents if s.startswith("ali") and s.endswith(".fits")]
  #imgs = [s for s in dir_contents if s.startswith("shftd") and s.endswith(".fits")]
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
  core = ((obs - 1.) / err)**2.
  return core.sum()/len(obs-1)

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

def nanmedian(arr, **kwargs):
    """
    Returns median ignoring NAN
    """
    return ma.median( ma.masked_where(arr!=arr, arr), **kwargs )

def bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,offset_l,offset_r):
	group = T_LC_R_C[group_size*i-offset_l:group_size*(i+1)-offset_r]
	binned.append(group.mean())
	binned_jd.append(JD_R[group_size*i-offset_l:group_size*(i+1)-offset_r].mean())
	uncertainties.append(	np.sqrt( 	(err_target_w_mean_norm[group_size*i-offset_l:group_size*(i+1)-offset_r]**2).sum())/(len(group)-1) 	)
	bin_err.append(group.std()/np.sqrt(len(group)))
	binned_err.append(np.array([uncertainties[i],bin_err[i]]).max())
	
	print group.mean(),"\t[",group_size*i-offset_l,":",group_size*(i+1)-offset_r,"] ",(group_size*(i+1)-offset_r)-(group_size*i-offset_l),"\twith uncertainty ",group.std()/np.sqrt(len(group)),"\tor ",uncertainties[i],"\tchosen: ",binned_err[i]

	f = open(folder+'/'+str(folder[23:])+'_bins.list', 'w+')
	print >> f, group.mean(),"\t[",group_size*i-offset_l,":",group_size*(i+1)-offset_r,"] ",(group_size*(i+1)-offset_r)-(group_size*i-offset_l),"\twith uncertainty ",group.std()/np.sqrt(len(group)),"\tor ",uncertainties[i],"\tchosen: ",binned_err[i]
	
def bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,offset_l,offset_r):
	group_REF = REF_LC_R_C[group_size*i-offset_l:group_size*(i+1)-offset_r]
	binned_REF.append(group_REF.mean())
	uncertainties_REF.append(	np.sqrt( 	(err_REF_w_mean_norm[group_size*i-offset_l:group_size*(i+1)-offset_r]**2).sum())/(len(group_REF)-1) 	)
	bin_err_REF.append(group_REF.std()/np.sqrt(len(group_REF)))
	binned_err_REF.append(np.array([uncertainties_REF[i],bin_err_REF[i]]).max())
 
	
	print group_REF.mean(),"\t[",group_size*i-offset_l,":",group_size*(i+1)-offset_r,"] ",(group_size*(i+1)-offset_r)-(group_size*i-offset_l),"\twith uncertainty ",group_REF.std()/np.sqrt(len(group_REF)),"\tor ",uncertainties_REF[i],"\tchosen: ",binned_err_REF[i]

	f = open(folder+'/'+str(folder[23:])+'_bins_REF.list', 'w+')
	print >> f, group_REF.mean(),"\t[",group_size*i-offset_l,":",group_size*(i+1)-offset_r,"] ",(group_size*(i+1)-offset_r)-(group_size*i-offset_l),"\twith uncertainty ",group_REF.std()/np.sqrt(len(group_REF)),"\tor ",uncertainties_REF[i],"\tchosen: ",binned_err_REF[i]
	
def do_plot(folder):
  ''' Read Julian Date and airmass info from .fits headers '''
  airmass, JD = get_airmass_JD(folder)
  airmass = np.array(airmass)
  JD = np.array(JD)
  fwhm = loadtxt(folder+'/data/fwhm.txt',usecols=(0,))
  
  ''' Calculate the seeing '''
  pixel_scale = 0.288#par.pixel_scale[0]
  seeing = fwhm*pixel_scale
  
  ''' Establishing directories '''
  dir_contents = listdir(folder)
  folders = [s for s in dir_contents if s.startswith("A")]
  data_folder = [s for s in dir_contents if s.startswith("data")]
  directory = folder+'/'+folders[0]
  data_directory = folder+'/'+data_folder[0]
  T = loadtxt(directory+'/xyf.list')
  
  ''' Convert JD days --> hours '''
  #JD_1 = JD[0]
  #JD = (JD-JD_1)*24
  
  ''' Defining the number of reference stars '''
  ds9 = loadtxt(folder+'/data/ds9.reg')
  global N_ref_stars
  N_ref_stars = len(ds9)-1
  
  ''' Creating 2D lists list[star][frame] '''
  flux_r = [[] for _ in range(N_ref_stars)]
  mag_r = [[] for _ in range(N_ref_stars)]
  err_r = [[] for _ in range(N_ref_stars)]
  err_r_MAG = [[] for _ in range(N_ref_stars)]
  flux_r_norm = [[] for _ in range(N_ref_stars)]
  err_r_norm = [[] for _ in range(N_ref_stars)]
  stddev = [[] for _ in range(N_ref_stars)]
  stddev_norm = [[] for _ in range(N_ref_stars)]
  
  ''' Obtaining the photometry data '''
  x, y, flux, mag, merr, msky, xcenter, ycenter = loadtxt(directory+'/xyf.list', unpack=True)
  x_target = x[::N_ref_stars+1]
  y_target = y[::N_ref_stars+1]
  
  ''' Reading the Target Info '''
  global flux_T
  flux_T = flux[::N_ref_stars+1]
  flux_T_MAG = mag[::N_ref_stars+1]
  
  ''' Normalising Target Flux '''
  flux_target_norm = flux_T/flux_T.mean()
  err_T = merr[::N_ref_stars+1]*flux_T/1.0857	# Convert from mag error to flux error
  err_T_MAG = merr[::N_ref_stars+1]+1e-6
  err_target_norm = err_T/flux_T.mean()
  
  ''' Reading reference star info and normalising '''
  for i in range(N_ref_stars):
	flux_r[i] = np.array(flux[i+1::N_ref_stars+1])
	mag_r[i] = np.array(mag[i+1::N_ref_stars+1])
	err_r[i] = np.array(merr[i+1::N_ref_stars+1]*flux_r[i]/1.0857)
	err_r_MAG[i] = np.array(merr[i+1::N_ref_stars+1] )
	''' Normalising the reference flux and errors '''
	flux_r_norm[i] = np.array(flux_r[i]/flux_r[i].mean())
	err_r_norm[i] = np.array(err_r[i]/flux_r[i].mean()+1e-12) #adding 1e-12 as IRAF has errors stopping at 0.0001 and I don't want to divide by zero.
	stddev[i] = flux_r[i].std()
	stddev_norm[i] = flux_r_norm[i].std()
  
  ''' Initialise the arrays for weighted mean calculation '''
  flux_r_norm = array(flux_r_norm)
  err_r_norm = array(err_r_norm)
  variance_norm = array(err_r_norm*err_r_norm)
  
  ''' Find the mean of the raw flux of the reference stars '''
  norm_ref_star_flux_sum = np.array([sum(a) for a in zip(*(flux_r_norm))])
  norm_ref_star_flux_sum = norm_ref_star_flux_sum/N_ref_stars
  
  fail_condition=1
  while (fail_condition==1):

    ''' Computing the weighted mean of the reference flux for each ref star'''
    norm_ref_star_flux_weight_sum=np.zeros(len(flux_r_norm[0]))
    sigma_weighted_mean_ref_stars=np.zeros(len(flux_r_norm[0]))
  
    for j in range(len(flux_r_norm[0])):
  	  wmean = 0.
	  total_weight = 0.
	  for i in range(N_ref_stars):  
	    wmean += (flux_r_norm[i][j] / variance_norm[i][j])
	    total_weight += (1. / variance_norm[i][j])

	  norm_ref_star_flux_weight_sum[j] = wmean / total_weight
	  sigma_weighted_mean_ref_stars[j] = np.sqrt(1./total_weight)


    print norm_ref_star_flux_weight_sum
    print sigma_weighted_mean_ref_stars
    sys.exit()
    
    #comparison ref star
    ref_number = 1
    
    
    ''' FINAL ERROR BARS '''
    err_target_w_mean_norm = np.sqrt(err_target_norm**2+sigma_weighted_mean_ref_stars**2)
    print err_target_w_mean_norm
    sys.exit()
    err_REF_w_mean_norm = np.sqrt(err_r_norm[ref_number]**2+sigma_weighted_mean_ref_stars**2)
    err_ref = np.sqrt(sigma_weighted_mean_ref_stars**2+err_r_norm**2)   # total error on the normalised and shape-corrected flux of each ref star
  
    ''' Printing info about target and ref stars to screen '''
    global f_ave,f_med,f_stddev, f_t_stddev,chi2_red,SN_ave,SN_med,f_t_norm,f_t_stddev
    f_ave = []
    f_med = []
    f_stddev = []
    chi2_red = []
    SN_ave = []
    SN_med = []

    for i in range(N_ref_stars):
	  f_med.append(np.median(flux_r[i]))
	  f_ave.append(flux_r[i].mean())
	  f_stddev.append((flux_r_norm[i]/norm_ref_star_flux_weight_sum).std())
	  chi2_red.append(chi2(flux_r_norm[i]/norm_ref_star_flux_weight_sum,err_ref[i]))
	  '''Need to use the error on the normalised flux to determine the S/N
      as mean Signal = 1 and noise is the error on that normalised flux value.
      SN_ave.append(1./err_r_MAG[i].mean())
	  SN_med.append(np.median(1./err_r_MAG[i]))'''
	  
	  SN_ave.append(1./err_r_norm[i].mean()) 
	  SN_med.append(np.median(1./err_r_norm[i]))
  
    SN_med, f_stddev = zip(*sorted(zip(SN_med,f_stddev))) #Sort by SN
    ''' Reject comparison stars based on difference between std. dev. of light-curve
    and median std. dev. of the individual points, and the reduced chi squared value 
    (assuming constant flux model),

    Reject by setting errors to very large value so that weighted mean gives them a 
    negligible weight. 

    If any star fails test, loop back and repeat calculation of weighted mean comp star. 
    '''

    fail_condition=0

    for i in range(N_ref_stars):
      print ''+str(i),'\t',int(np.median(flux_r[i])),'\t',int(flux_r[i].mean()),'\t',round((flux_r_norm[i]/norm_ref_star_flux_weight_sum).std(),4),'\t',round(chi2(flux_r_norm[i]/norm_ref_star_flux_weight_sum,err_ref[i]),2),'\t',round(1/err_r_MAG[i].mean(),2),'\t\t',round(np.median(1/err_r_MAG[i]),2)    #,'\t',round(err_r_MAG[i].mean(),4)

    for i in range(N_ref_stars):
      if (chi2_red[i] > 1000.0) or ((f_stddev[i]/(1./SN_med[i]) ) > 1000.00) : # 4 then three
	if not (err_r_norm[i][0]==1000.): 
          err_r_norm[i]=1000. 
	  variance=err_r_norm * err_r_norm
          fail_condition=1
        print i, 'bad ref star...'

    '''if stats are okay then exit loop '''
    if fail_condition==0:
      print 'breaking'
      break

  # Target light curve:
  
  T_LC = flux_target_norm/norm_ref_star_flux_weight_sum
  
  # Sample ref light curve
  
  REF_LC = flux_r_norm[ref_number]/norm_ref_star_flux_weight_sum


  # Removing outliers
  MAD_LC = MAD(T_LC, c=0.6745, axis=None)
  
  removed_jd_points = []
  removed_flux_points = []
  removed_fwhm_points = []
  
  print len(T_LC),len(fwhm)
  i_removed = []
  cut_off = 30.0 # in multiples of MAD
  for i in range(len(T_LC)):
    if (T_LC[i] < (1.0-cut_off*MAD_LC)) or (T_LC[i] >= (1.0+cut_off*MAD_LC)):
      print JD[i],"-->\t",T_LC[i],"\t",fwhm[i]
      removed_jd_points.append(JD[i])
      removed_flux_points.append(T_LC[i])
      removed_fwhm_points.append(fwhm[i])
      i_removed.append(i)

  print "Removed ",len(removed_jd_points)," points.\n"
  print i_removed


  # END removing outliers
  T_LC_R = array([x for x in list(T_LC) if x not in removed_flux_points])
  REF_LC_R = np.delete(REF_LC, i_removed)
  JD_R = array([x for x in list(JD) if x not in removed_jd_points])
  fwhm_R = [x for x in list(fwhm) if x not in removed_fwhm_points]
  '''
  fwhm_R.append(5.972)

  fwhm_R.append(5.972)
  fwhm_R.append(5.972)
  fwhm_R.append(5.972)
  fwhm_R.append(5.972)
  fwhm_R.append(3.294)

  #fwhm_R.append(5.972)
  '''


  fwhm_R = array(fwhm_R)
  
  print len(JD_R),len(fwhm_R),len(T_LC_R)
  slope, intercept, r_value, p_value, std_err = stats.linregress(fwhm_R,T_LC_R)
  
  line = slope*fwhm_R+intercept
  seeing = fwhm_R*pixel_scale

  print 'Standard Deviation', std_err

  print "\n Star:\tF_med:\tF_ave:\tSTDEV:\tChi^2:\tS/N_ave:\tS/N_med:"     #MERR:\n"
  print " T","\t",int(np.median(flux_T)),"\t",int(flux_T.mean()),'\t',round((T_LC).std(),3),'\t',round(chi2(T_LC,err_target_w_mean_norm),1),'\t',round(1/err_T_MAG.mean(),1), '\t\t', round(np.median(1/err_T_MAG),1)    #, '\t',round((err_target_w_mean_norm*1.0857/flux_target_norm).mean(),4)
    

  # Creating a list of all reference stars (see ref_stars.list)
  f = open(folder+'/'+str(folder[23:])+'_ref_stars.list', 'w+')
  print >> f, "Star:\tF_med:\tF_ave:\tSTDEV:\tChi^2:\tS/N_ave:\tS/N_med:"     #MERR:\n"
  for i in range(N_ref_stars):
	print ' '+str(i),'\t',int(np.median(flux_r[i])),'\t',int(flux_r[i].mean()),'\t',round((flux_r_norm[i]/norm_ref_star_flux_weight_sum).std(),4),'\t',round(chi2(flux_r_norm[i]/norm_ref_star_flux_weight_sum,err_ref[i]),2),'\t',round(1/err_r_MAG[i].mean(),2),'\t\t',round(np.median(1/err_r_MAG[i]),2)    #,'\t',round(err_r_MAG[i].mean(),4)
	print >> f,i,'\t',int(np.median(flux_r[i])),'\t',int(flux_r[i].mean()),'\t',round((flux_r_norm[i]/norm_ref_star_flux_weight_sum).std(),4),'\t',round(chi2(flux_r_norm[i]/norm_ref_star_flux_weight_sum,err_ref[i]),2),'\t',round(1/err_r_MAG[i].mean(),2),'\t\t',round(np.median(1/err_r_MAG[i]),2)     #,'\t',1/err_r_MAG[i].mean()
  f.close()
  print ""


  #f = open('stats.txt', 'a')
  #print >> f,str(folder[28:]),len(flux_target_norm),"\t",MAD(T_LC_R, c=0.6745, axis=None),"\t",flux_T.mean(),"\t",err_T.mean()
  #f.close()
  
  # Step C
  global binned, bin_err, binned_err, binned_jd, uncertainties
  global binned_REF, bin_err_REF, binned_err_REF, uncertainties_REF
  

  binned = []
  bin_err = []
  binned_err = []
  binned_jd = []
  uncertainties = []
  
  binned_REF = []
  bin_err_REF = []
  binned_err_REF = []
  uncertainties_REF = []
  
  group_size = 20                                                       # Change this
  T_LC_R_C = T_LC_R/line     # Detrending taken into account
  REF_LC_R_C = REF_LC_R/line
  print "\t\t\t\t\t\t  sig/sqrt(n)\t\t\terr"


  for i in range(8):                                                  # Change this too
    #print i,"\n"
    #bin_points(i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,0,0)
    #''' 
    if i == [0,1,2,3,4,5]:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,0,0)
    else:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,0,7)
    
    '''
    elif i == 1:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,4,5)
    elif i == 2:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,5,5)
    elif i == 3:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,5,5)
    elif i == 4:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,5,5)
    elif i == 5:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,5,5)
    elif i == 6:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,5,5)
    elif i == 7:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,5,6)
    elif i == 8:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,6,6)
    elif i == 9:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,6,6)
    elif i == 10:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,6,12)
    elif i == 11:
      bin_points(folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,12,-8)
    '''
  
  print "\n\nREF\n"    
  for i in range(8):                                                  # Change this too
    #print i,"\n"
    #bin_points(i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,0,0)
    #folder,i,group_size,err_target_w_mean_norm,JD_R,T_LC_R_C,offset_l,offset_r
    #''' 
    if i == [0,1,2,3,4,5]:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,0,0)
    else:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,0,7)
    
    '''
    elif i == 1:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,4,5)
    elif i == 2:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,5,5)
    elif i == 3:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,5,5)
    elif i == 4:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,5,5)
    elif i == 5:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,5,5)
    elif i == 6:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,5,5)
    elif i == 7:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,5,6)
    elif i == 8:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,6,6)
    elif i == 9:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,6,6)
    elif i == 10:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,6,12)
    elif i == 11:
      bin_points_REF(folder,i,group_size,err_REF_w_mean_norm,JD_R,REF_LC_R_C,12,-8)
    '''

  f.close()

  print "\nNumber of Light Curve points: ",len(T_LC)

  print "\nChi^2_red:\t",chi2(np.array(binned),np.array(binned_err)),"\tDOF:\t",len(binned)-1
  print "Due to chance [%] ",(1-stats.chi2.cdf(chi2(np.array(binned),np.array(binned_err))*(len(binned)-1), len(binned)-1))*100.
  
  print "\n",chi2(np.array(binned),np.array(binned_err)),"\t",len(binned)-1,"\t",(1-stats.chi2.cdf(chi2(np.array(binned),np.array(binned_err))*(len(binned)-1), len(binned)-1))*100.,"\t",np.mean(bin_err),"\t",np.mean(uncertainties),"\t",(norm_ref_star_flux_sum/norm_ref_star_flux_weight_sum).std()

  # =============================  PLOTS  =============================
  plt.figure(1)
  ax = subplot('211')
  title('Star: '+str(folder[23:]))
  plot([0,JD_R[-1]],[1,1],'--r')
  errorbar(binned_jd,binned,yerr=binned_err,fmt='.',elinewidth=1.0,capsize=3,markersize=5,color="black")
  text(0.1, 1.01, r"$\chi^2_{\rm{red}}=$"+str(round(chi2(np.array(binned),np.array(binned_err)),1)), fontsize=16, color='black')
  
  f = open('all_data/'+folder[23:]+'/target_lc/'+folder[23:]+'_LC_binned.txt', 'w+')
  for j in range(len(binned_jd)):
    print >> f, binned_jd[j],binned[j],binned_err[j]
  f.close() 

  f = open('all_data/'+folder[23:]+'/refs/'+folder[23:]+'REF_bin.txt', 'w+')
  for j in range(len(binned_jd)):
    print >> f, binned_jd[j],binned_REF[j],binned_err_REF[j]
  f.close() 
  
  xlabel('JD [hours]')
  ylabel('Relative Flux')
  #ylim(0.97,1.03)
  ax2 = subplot('212')
  plot(fwhm_R,line,'r-',fwhm_R,T_LC_R,'o')
  xlabel('FWHM')
  ylabel('Relative Flux')
  savefig(folder[23:]+'_trend.pdf')

  plt.figure(2)
  ax = subplot('311')
  title('Star: '+str(folder[23:]))
  #errorbar(JD,flux_target_norm/norm_ref_star_flux_weight_sum, yerr=err_target_w_mean_norm, fmt='.',elinewidth=1.0,capsize=3,markersize=5,color="black")
  plot(JD_R,T_LC_R,'.k',alpha=0.2)
  #plot(JD_R,line,'-r')
  plot(JD_R,T_LC_R_C,'.r')
  #plot(JD,flux_target_norm/flux_r_norm[1],'.k')
  
  f = open('all_data/'+folder[23:]+'/target_lc/'+folder[23:]+'_LC_.txt', 'w+')
  for j in range(len(T_LC_R_C)):
    print >> f, JD_R[j],T_LC_R_C[j],'\t',err_target_w_mean_norm[j]
  f.close()  

  ylabel('Relative Flux')
  #ylim(0.96,1.04)
  
  ax2 = subplot('312')
  plot(JD,np.ones(len(JD)),'-',color="black",alpha=1.0)
  
  f = open('all_data/'+folder[23:]+'/refs/'+folder[23:]+'_REF.txt', 'w+')
  for j in range(len(T_LC_R_C)):
    print >> f, JD_R[j],REF_LC_R_C[j],'\t',err_REF_w_mean_norm[j]
  f.close()
  
  for i in range(N_ref_stars):
  	
  	'''
  	# Save ref star light curves to a .txt file
  	f = open('all_data/'+folder[23:]+'/refs/'+folder[23:]+'_ref_'+str(i)+'.txt', 'w+')
	for j in range(len(JD)):
	  #print >> f, JD[j], flux_r_norm[i][j]/norm_ref_star_flux_weight_sum[j],'\t',err_ref[i][j]
	  print >> f, JD_R[j], REF_LC_R_C[j],'\t',err_REF_w_mean_norm[j]
	f.close()	
	'''
	
	# Plot the ref star light curves
	plot(JD,flux_r_norm[i]/norm_ref_star_flux_weight_sum,'.',alpha=0.3)

  #plot(JD,flux_r_norm[0]/norm_ref_star_flux_weight_sum,'.',alpha=0.3)
  errorbar(JD,norm_ref_star_flux_sum/norm_ref_star_flux_weight_sum,yerr=sigma_weighted_mean_ref_stars, fmt='.',elinewidth=1.0,capsize=3,markersize=5,color="black")
  #ylim(0.98,1.02)
  
  # Uncomment below if you want to look at a particular reference star
  #errorbar(JD,flux_r_norm[0]/norm_ref_star_flux_weight_sum, yerr=err_ref[0], fmt='.',elinewidth=1.0,capsize=3,markersize=5,color="red")
  
  #legend(('0', '1', '2', '3', '4', '5','6','7','8','9','10', '11', '12', '13', '14', '15','16','17','18','19','20','21'),'upper left', shadow=True)
  ylabel('Reference stars')
  ylim(0.90,1.10)
  
  ax3 = subplot('313')
  plot(JD_R,seeing,'.b',alpha=0.5)
  plot(JD,norm_ref_star_flux_weight_sum,'.r',alpha=0.5)
  plot(JD,airmass,'--g',alpha=1.0)
  ylabel('Seeing (Blue) / Flux (Red)')
  xlabel('Hours')
  #print data_directory+'/'+imgs[5]
  #header = pyfits.open(data_directory+'/'+imgs[5])
  #print header[0].header['OBJECT']
  #print folder[26:]
  #print "Creating: ",folder[26:]+'.pdf'
  savefig(folder[23:]+'.pdf')#,dpi=300)
  #ax.cla()
  #ax2.cla()
  #ax3.cla()

  # =============================  PLOTS END  =============================

  #show()
  
  print T_LC_R_C.std()

