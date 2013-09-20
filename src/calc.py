import numpy as np
from itertools import izip_longest
from pylab import *
import math
from scipy import stats

import unix as u
import statistical as st
import param as par

import datetime
import os

def plot_them(iteration,obj_name,JD,T_LC,err_target_w_mean_norm,REF_LCs,master_ref,err_r_norm_sum,airmass,np_airmass,poly2d_AIRMASS,N_ref_stars,err_ref_Rs,robust,fwhm):
 
  rainbow_colors = iter(cm.rainbow(np.linspace(0, 1, N_ref_stars)))
  
  if (iteration == 0):
    target_colour = '#D8D8D8'
    master_ref_colour = '#D8D8D8'
  else:
    target_colour = 'red'
    master_ref_colour = 'red'
 
  fig=plt.figure(1,figsize=(8.27,11.69))
  
  # Light curve plot
  main_panel=fig.add_subplot(3,2,1)
  main_panel.cla()
  plt.errorbar(JD,T_LC,yerr=err_target_w_mean_norm,fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color=target_colour, ecolor="black")
  ylabel('Relative Flux')
  
  # Master ref plot
  middle_panel=fig.add_subplot(3,2,3,sharex=main_panel, sharey=main_panel)
  middle_panel.cla()
  for i in range(len(REF_LCs)):
    plt.plot(JD,REF_LCs[i],'.k',alpha=0.1)
  plt.errorbar(JD,master_ref,yerr=err_r_norm_sum,fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color=master_ref_colour, ecolor="black")
  ylabel('Relative Flux')
  plt.ylim(T_LC.min()-0.01,T_LC.max()+0.01)
  
  # Airmass trend plot
  if (iteration == 0):
    bottom_panel=fig.add_subplot(3,2,5)
    bottom_panel.cla()
    plt.plot(airmass,master_ref,'ok')
    plt.plot(np_airmass,poly2d_AIRMASS,'-r')
    xlabel('Airmass')
    ylabel('Relative Flux')
  else:
    bottom_panel=fig.add_subplot(3,2,5)
    bottom_panel.cla()
    #plt.plot(airmass,master_ref,'ok')
    plt.plot(JD,fwhm,'ok')
    xlabel('Days')
    ylabel('FWHM')
    
  # Ref star plots on the right hand side
  
  if (iteration == 0):
    for i in range(N_ref_stars):
      side_panel = fig.add_subplot(N_ref_stars,2,2*(i+1),sharex=main_panel, sharey=main_panel)
      plt.errorbar(JD, REF_LCs[i],yerr=err_ref_Rs[i],fmt='.',elinewidth=1.3,capsize=3,markersize=12,ecolor='black',color=next(rainbow_colors),markeredgecolor='black')

  else:
    for i in range(N_ref_stars):
      side_panel = fig.add_subplot(N_ref_stars,2,2*(i+1),sharex=main_panel, sharey=main_panel)
      side_panel.cla() 
      plt.errorbar(JD, REF_LCs[i],yerr=err_ref_Rs[i],fmt='.',elinewidth=1.3,capsize=3,markersize=12,ecolor='black',color=next(rainbow_colors),markeredgecolor='black')
      plt.ylim(np.array([T_LC.min(),master_ref.min()]).min()-0.01,np.array([T_LC.max(),master_ref.max()]).max()+0.01)
  
  
  if par.variable_aperture[0] == "yes":
    aperture = "variable"
  else:
    aperture = "constant"
  minorticks_on()
  if robust > 1.0 and iteration > 0:
    plt.savefig('plots/'+'V_'+str(iteration)+'_'+obj_name+'_'+aperture+'_aperture.pdf')
  else:
    plt.savefig('plots/'+str(iteration)+'_'+obj_name+'_'+aperture+'_aperture.pdf')
  plt.draw()
  if (iteration == 0):
    clf()

def calc(obj_name,N_ref_stars,N_bad_refs,airmass,np_airmass,poly2d_AIRMASS,fwhm,JD,T_LC,err_target_w_mean_norm,norm_ref_star_flux_weight_sum_mean,REF_LCs,err_r_norm_sum,err_ref_Rs,master_ref,iteration,N_orig_stars,ds9_file):
  
  MAD_factor = par.cut_offs[1]
  MAD_cutoff = []
  STD_cutoff = []
   
  for i in range(len(REF_LCs)):
    MAD_cutoff.append(st.MAD(REF_LCs[i]))
    STD_cutoff.append(REF_LCs[i].std())
  
  MAD_median = np.median(MAD_cutoff)
  STD_cutoff = np.median(STD_cutoff)
  MAD_STD = np.std(MAD_cutoff)
  
  ref_std = []
  for i in range(len(REF_LCs)):
    ref_std.append(REF_LCs[i].std())
  ref_std = np.array(ref_std)
  mean_ref_std = ref_std.mean()
  
  robust = st.robust(T_LC,err_target_w_mean_norm)
  chi2_red = st.chi2(T_LC,err_target_w_mean_norm)
  chi2_master_ref = st.chi2_master_ref(T_LC,err_target_w_mean_norm,master_ref,err_r_norm_sum)
  for i in range(len(T_LC)):
    if T_LC[i] == T_LC.max():
      top_err = err_target_w_mean_norm[i]
    if T_LC[i] == T_LC.min():
      bottom_err = err_target_w_mean_norm[i]
  #amplitude = (T_LC.max()-T_LC.min())*100.," +- ",round(np.sqrt(top_err**2+bottom_err**2)*100.,1
  #amplitude_print = round((T_LC.max()-T_LC.min())*100.,2)," +- ",round(np.sqrt(top_err**2+bottom_err**2)*100.,1)
  amplitude = (T_LC.max()-T_LC.min())*100.
  amplitude_err = np.sqrt(top_err**2+bottom_err**2)*100.
  
  # Print to screen
  print "\n"
  print "\tChi2:\t",round(chi2_master_ref,2)
  print "\tRobust:\t",round(robust,2)
  print "\tMean Ref std:\t\t",round(mean_ref_std,5)
  print "\tSTDDEV Target:\t\t",round(master_ref.std(),5)
  print "\tSTDDEV Master Ref:\t",round(T_LC.std(),5)
  print "\tAmplitude:\t\t",round(amplitude,2),"+-",round(amplitude_err,2),"\n"
  
  if (iteration == 0):
    print "\tRef Stars:\t\t","cut-off: "+str(round(MAD_median+MAD_factor*MAD_STD,5))+"\n"
    print "\tChi2:","\t\t","Robust:","\t","MAD:","\t\t","ADU_max"
    print "\t---------------------------------------------------"
    for i in range(len(REF_LCs)):
      #print "\t",round(st.chi2(REF_LCs[i],err_ref_Rs[i]),2),'\t\t',round(st.robust(REF_LCs[i],err_ref_Rs[i]),2),'\t\t',round(st.MAD(REF_LCs[i]),5)


      if ((st.MAD(REF_LCs[i]) > MAD_median+MAD_factor*MAD_STD) or (st.chi2(REF_LCs[i],err_ref_Rs[i]) > 50.) or (ds9_file[4][i+1]>par.cut_offs[5]) or (ds9_file[4][i] < par.cut_offs[4])):
        print "\t",round(st.chi2(REF_LCs[i],err_ref_Rs[i]),2),'\t\t',round(st.robust(REF_LCs[i],err_ref_Rs[i]),2),'\t\t',round(st.MAD(REF_LCs[i]),5),"\t",ds9_file[4][i+1],u.red('  Ref: '+str(i))
      else:
        print "\t",round(st.chi2(REF_LCs[i],err_ref_Rs[i]),2),'\t\t',round(st.robust(REF_LCs[i],err_ref_Rs[i]),2),'\t\t',round(st.MAD(REF_LCs[i]),5),"\t",ds9_file[4][i+1]   

    f = file("temp.txt", "w+")    
    for i in range(len(REF_LCs)):
      if ((st.MAD(REF_LCs[i]) > MAD_median+MAD_factor*MAD_STD) or (st.chi2(REF_LCs[i],err_ref_Rs[i]) > 50.)):
        print >> f, i  
    f.close()
    print "\t----------------------------------------" 

  
  if (iteration == 1):
    DOF = len(master_ref)-1
    chance_variable = round((1-stats.chi2.cdf(chi2_master_ref*DOF, DOF))*100.,4)
    if (T_LC.std() > master_ref.std() and robust >= 1.0 and chi2_master_ref >= 2.0 and chance_variable <= 5.0):
      variable = 'y'
    else:
      variable = 'n'


    f = open('data.txt', 'a+')
    print >> f, obj_name,'\t',round(chi2_red,2),'\t',round(chi2_master_ref,2),'\t',round(robust,2),'\t',round(T_LC.std(),6),'\t',round(mean_ref_std,6),'\t',round(master_ref.std(),6),'\t',str(N_orig_stars),'\t',str(N_bad_refs),'\t',str(N_ref_stars),'\t',DOF,'\t',round(amplitude,2),'\t',round(amplitude_err,2),'\t',chance_variable,'\t',variable
    f.close()
    
    f = open('plots/T-Y_Plots/LCs/data/W0458.txt', 'w+')
    for i in range(len(T_LC)):
      print >> f, JD[i],T_LC[i],err_target_w_mean_norm[i],master_ref[i],err_r_norm_sum[i]
    
    f.close()

    #Calculating Photometric Quality
    #phot_qual = []
    #for i in range(len(T_LC)):
    #  phot_qual.append(np.median(err_target_w_mean_norm[i]))
    
    #f = open('../plots/L-T_Plots/Histogram_phot/data/phot_qual_binned.txt', 'a+')
    #print >> f, np.median(err_target_w_mean_norm)  
    #f.close()
  
  #obj_name,JD,T_LC,err_target_w_mean_norm,REF_LCs,master_ref,err_r_norm_sum,airmass,np_airmass,poly2d_AIRMASS,N_ref_stars,err_ref_Rs,robust,fwhm
  plot_them(iteration,obj_name,JD,T_LC,err_target_w_mean_norm,REF_LCs,master_ref,err_r_norm_sum,airmass,np_airmass,poly2d_AIRMASS,N_ref_stars,err_ref_Rs,robust,fwhm)

