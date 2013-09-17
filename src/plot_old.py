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
import param as par
from os import listdir
import time
import matplotlib
matplotlib.rcParams.update({'font.size': 9})

import xml.etree.cElementTree as ET

from matplotlib.dates import date2num
from datetime import datetime

def get_airmass_JD():
  dir_contents = listdir('data')
  imgs = [s for s in dir_contents if s.startswith("ali") and s.endswith(".fits")]
  imgs = sorted(imgs)
  headers = []
  airmass = []
  JD = []
  for i in range(len(imgs)):
    headers.append(pyfits.open('data/'+imgs[i]))
    JD.append(headers[i][0].header['MJD-OBS'])
    airmass.append(headers[i][0].header['HIERARCH ESO TEL AIRM START'])
  return airmass, JD

def chi2(obs,err):
  core = ((obs - 1.) / err)**2.
  return core.sum()/len(obs-1)


  '''
  target = ET.SubElement(root, obj)
  flux = ET.SubElement(target, "flux")
  average = ET.SubElement(flux, "average", attrib={"hello"})
  average.text = str(f_ave[i])
  median = ET.SubElement(flux, "median")
  median.text = str(f_med[i])
  stats = ET.SubElement(target, "stats")
  stddev = ET.SubElement(stats, "stddev")
  stddev.text = str(f_stddev[i])
  chi2 = ET.SubElement(stats, "chi2")
  chi2.text = str(chi2_red[i])
  SN = ET.SubElement(target, "SN")
  SN_a = ET.SubElement(SN, "average")
  SN_a.text = str(SN_ave[i])
  SN_m = ET.SubElement(SN, "median")
  SN_m.text = str(SN_med[i])
  '''


def output_data():
  global dataset_mean, dataset_med, categories
  dataset_mean = ET.SubElement(parent, 'dataset', seriesName="Mean", color="346BAB")
  dataset_med = ET.SubElement(parent, 'dataset', seriesName="Median", color="B23732")
  categories = ET.SubElement(parent, 'categories')
  for i in range(N_ref_stars):
    child("Ref"+str(i),i)
  target = ET.SubElement(parent, 'trendlines')
  ET.SubElement(target, 'line', startValue=str(np.median(flux_T)),displayValue='Target Star', thickness='3',color='E00000')
  tree = ET.ElementTree(parent)
  tree.write("src/xml/data.xml")
  return

def lc_data():
  global dataset_target, dataset_ref, categories
  dataset_target = ET.SubElement(parent_lc, 'dataset', seriesName="Target", color="346BAB")
  dataset_ref = ET.SubElement(parent_lc, 'dataset', seriesName="Median", color="346BAB")
  categories = ET.SubElement(parent_lc, 'categories')
  for i in range(len(JD)):
    child_lc("Ref"+str(i),i)
  target = ET.SubElement(parent_lc, 'trendlines')
  ET.SubElement(target, 'line', startValue=str(1.), thickness='3',color='E00000',showOnTop='0')
  tree = ET.ElementTree(parent_lc)
  tree.write("src/xml/lc.xml")
  return

def child(obj,i):
  ET.SubElement(dataset_mean, 'set', name=obj, value=str(f_ave[i]), color="346BAB", FontSize="19")
  ET.SubElement(dataset_med, 'set', name=obj, value=str(f_med[i]), color="B23732")
  ET.SubElement(categories, 'category', name=obj)

def child_lc(obj,i):
  ET.SubElement(dataset_target, 'set', name=str(JD[i]), value=str(flux_target_norm[i]/norm_ref_star_flux_weight_sum[i]), color="346BAB", FontSize="19",alpha='0.01')
  #ET.SubElement(dataset_mean, 'set', name=obj, value=str(f_med[i]), color="B23732")
  ET.SubElement(categories, 'category', name=obj)

''' Read Julian Date and airmass info from .fits headers '''
airmass, JD = get_airmass_JD()
airmass = np.array(airmass)
JD = np.array(JD)
fwhm = loadtxt('data/fwhm.txt',usecols=(0,))

''' Calculate the seeing '''
pixel_scale = par.pixel_scale[0]
seeing = fwhm*pixel_scale

''' Defines which folder you are plotting. See param.py '''
name  = par.object_name[0]
folder = par.folder[0]
T = loadtxt(folder+'/xyf.list')

''' Convert JD days --> hours '''
JD_1 = JD[0]
JD = (JD-JD_1)*24

''' Defining the number of reference stars '''
ds9 = loadtxt('ds9.reg')
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
x, y, flux, mag, merr = loadtxt(folder+'/xyf.list', unpack=True)
x_target = x[::N_ref_stars+1]
y_target = y[::N_ref_stars+1]

''' Reading the Target Info '''
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
  err_r_norm[i] = np.array(err_r[i]/flux_r[i].mean()+1e-6) #adding 1e-6 as IRAF has errors stopping at 0.0001 and I don't want to divide by zero.
  stddev[i] = flux_r[i].std()
  stddev_norm[i] = flux_r_norm[i].std()

''' Initialise the arrays for weighted mean calculation '''
flux_r_norm = array(flux_r_norm)
err_r_norm = array(err_r_norm)
variance_norm = array(err_r_norm*err_r_norm)

''' Find the mean of the raw flux of the reference stars '''
norm_ref_star_flux_sum = np.array([sum(a) for a in zip(*(flux_r_norm))])
norm_ref_star_flux_sum = norm_ref_star_flux_sum/N_ref_stars

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

norm_ref_star_flux_weight_sum = array(norm_ref_star_flux_weight_sum)
sigma_weighted_mean_ref_stars = array(sigma_weighted_mean_ref_stars)


''' FINAL ERROR BARS '''
err_ref = np.sqrt(sigma_weighted_mean_ref_stars**2+err_r_norm**2)   # total error on the normalised and shape-corrected flux of each ref star
err_target_w_mean_norm = np.sqrt(err_target_norm**2+sigma_weighted_mean_ref_stars**2)

''' Printing info about target and ref stars to screen '''

f_ave = []
f_med = []
f_stddev = []
chi2_red = []
SN_ave = []
SN_med = []
for i in range(N_ref_stars):
  f_ave.append(np.median(flux_r[i]))
  f_med.append(flux_r[i].mean())
  f_stddev.append((flux_r_norm[i]/norm_ref_star_flux_weight_sum).std())
  chi2_red.append(chi2(flux_r_norm[i]/norm_ref_star_flux_weight_sum,err_ref[i]))
  SN_ave.append(1./err_r_MAG[i].mean())
  SN_med.append(np.median(1./err_r_MAG[i]))

'''
print "\n Star:\tF_med:\tF_ave:\tSTDEV:\tChi^2:\tS/N_ave:\tS/N_med:"     #MERR:\n"
print " T(N)","\t",int(np.median(flux_T)),"\t",int(flux_T.mean()),'\t',round((flux_target_norm/norm_ref_star_flux_weight_sum).std(),3),'\t',round(chi2(flux_target_norm/norm_ref_star_flux_weight_sum,err_target_w_mean_norm),1),'\t',round(1/err_T_MAG.mean(),1), '\t\t', round(np.median(1/err_T_MAG),1)    #, '\t',round((err_target_w_mean_norm*1.0857/flux_target_norm).mean(),4)

round(1/err_T_MAG.mean(),1)
round(1/err_r_MAG[i].mean(),2)


# Creating a list of all reference stars (see ref_stars.list)
f = open('ref_stars.list', 'w+')
print >> f, "Star:\tF_med:\tF_ave:\tSTDEV:\tChi^2:\tS/N_ave:\tS/N_med:"     #MERR:\n"
for i in range(N_ref_stars):
  print ' '+str(i),'\t',int(np.median(flux_r[i])),'\t',int(flux_r[i].mean()),'\t',round((flux_r_norm[i]/norm_ref_star_flux_weight_sum).std(),4),'\t',round(chi2(flux_r_norm[i]/norm_ref_star_flux_weight_sum,err_ref[i]),2),'\t',round(1/err_r_MAG[i].mean(),2),'\t\t',round(np.median(1/err_r_MAG[i]),2)    #,'\t',round(err_r_MAG[i].mean(),4)
  print >> f,i,'\t',int(np.median(flux_r[i])),'\t',int(flux_r[i].mean()),'\t',round((flux_r_norm[i]/norm_ref_star_flux_weight_sum).std(),4),'\t',round(chi2(flux_r_norm[i]/norm_ref_star_flux_weight_sum,err_ref[i]),2),'\t',round(1/err_r_MAG[i].mean(),2),'\t\t',round(np.median(1/err_r_MAG[i]),2)     #,'\t',1/err_r_MAG[i].mean()
f.close()
'''

''' Plotting the results '''

''' Top Plot '''
subplot('311')
title('OB15-Red')
errorbar(JD,flux_target_norm/norm_ref_star_flux_weight_sum, yerr=err_target_w_mean_norm, fmt='.',elinewidth=1.0,capsize=3,markersize=5,color="black")
ylabel('Relative Flux')
#ylim(0.98,1.02)

''' Middle Plot '''
subplot('312')
plot(JD,np.ones(len(JD)),'-',color="black",alpha=1.0)

for i in range(N_ref_stars):
  plot(JD,flux_r_norm[i]/norm_ref_star_flux_weight_sum,'.',alpha=0.3)


errorbar(JD,norm_ref_star_flux_sum/norm_ref_star_flux_weight_sum,yerr=sigma_weighted_mean_ref_stars, fmt='.',elinewidth=1.0,capsize=3,markersize=5,color="black")

# Uncomment below if you want to look at a particular reference star
#errorbar(JD,flux_r_norm[0]/norm_ref_star_flux_weight_sum, yerr=err_ref[0], fmt='.',elinewidth=1.0,capsize=3,markersize=5,color="red")

#legend(('0', '1', '2', '3', '4', '5','6','7','8','9','10', '11', '12', '13', '14', '15','16','17','18','19','20','21'),'upper left', shadow=True)
ylabel('Reference stars')
#ylim(0.99,1.01)

''' Bottom Plot '''
subplot('313')
plot(JD,seeing,'.b',alpha=0.5)
plot(JD,norm_ref_star_flux_weight_sum,'.r',alpha=0.5)
plot(JD,airmass,'--g',alpha=1.0)
ylabel('Seeing (Blue) / Flux (Red)')
xlabel('Hours')
#savefig('OBJ48.pdf')#,dpi=300)
#show()


''' Creating Systematic Plot '''
'''
fig = figure()
fig.subplots_adjust(wspace = 0.35,hspace = 0.45)
fig.suptitle('Systematics')
ET
ax = fig.add_subplot(211)
ax.plot(flux_T.mean(),err_T.mean(),'*',markersize=15,color='black',markeredgecolor='black',markeredgewidth=1.0,alpha=1.0)
#ax.plot(np.arange(35000.),np.sqrt(np.arange(35000.)),'-r')
#ax.plot(np.arange(35000.),np.ones(35000.)*det_limit,'-r')
for i in range(N_ref_stars):
  ax.plot(flux_r[i].mean(),err_r[i].mean(),'.',markersize=8,color='red',markeredgecolor='black',alpha=0.9)
  #ax.plot(1/err_r_MAG[i].mean(),err_r[i].mean(),'.',markersize=8,color='red',markeredgecolor='black',alpha=0.9)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('$\\rm{Flux}$',fontsize=18)
ax.set_ylabel('$\sigma_{\\rm{Flux}}$',fontsize=18)
#savefig('Ross458C-systematics.pdf')


ax = fig.add_subplot(212)
ax.plot(1/err_T_MAG.mean(),err_T.mean(),'*',markersize=15,color='black',markeredgecolor='black',markeredgewidth=1.0,alpha=1.0)
for i in range(N_ref_stars):
  #ax.plot(flux_r[i].mean(),err_r[i].mean(),'.',markersize=8,color='red',markeredgecolor='black',alpha=0.9)
  ax.plot(1/err_r_MAG[i].mean(),err_r[i].mean(),'.',markersize=8,color='red',markeredgecolor='black',alpha=0.9)

 
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_xlabel('$\\rm{S/N}$',fontsize=18)
ax.set_ylabel('$\sigma_{\\rm{Flux}}$',fontsize=18)
'''
#show()
parent = ET.Element('graph')#,yAxisName='Revenue' )
parent.set('yAxisName', 'Flux')
parent.set('baseFontSize', '12')

parent_lc = ET.Element('graph')#,yAxisName='Revenue' )
parent_lc.set('yAxisName', 'Flux')
parent_lc.set('baseFontSize', '12')
parent_lc.set('showValues', '0')
parent_lc.set('shownames', '0')
parent_lc.set('yAxisMinValue', '0.85')
parent_lc.set('yAxisMaxValue', '1.15')
parent_lc.set('showShadow', '0')
parent_lc.set('anchorRadius', '3')
parent_lc.set('anchorSides','20')
parent_lc.set('anchorBorderThickness','2')
parent_lc.set('anchorBorderColor ','0000A0')
parent_lc.set('anchorBgColor','ADD8E6')
output_data()
lc_data()

