# To make the two panel plot used for the unbinned data

  '''
  if (iteration != 0):
    fig=plt.figure(2,figsize=(11.69,8.27))
    main_panel=fig.add_subplot(2,1,1)
    plt.errorbar(JD,T_LC,yerr=err_target_w_mean_norm,fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color=target_colour, ecolor="black")
    
    xlabel('Time [Hours]', fontsize=14)
    ylabel('Normalised Ref. Flux', fontsize=14)
    plt.ylim(0.97,1.02)
    
    bottom_panel=fig.add_subplot(2,1,2)
    plt.errorbar(JD,master_ref,yerr=err_r_norm_sum,fmt='.',elinewidth=1.3,capsize=3,markersize=12,markeredgecolor='black',color=target_colour, ecolor="black")
    xlabel('Time [Hours]', fontsize=14)
    ylabel('Normalised Ref. Flux', fontsize=14)
    plt.ylim(0.97,1.02)
  
  '''
  
	  '''
	  if (iteration == 0):
	    for i in range(N_ref_stars):
	      weights_Rs = weights
	      weighted_mean_Rs = weighted_mean
	      
	      weights_Rs = [x for x in weights_Rs if x not in weights_Rs[i]] # Removing the contribution by the comparison star
	      #weights_Rs_sum = np.array([np.median(a) for a in zip(*(weights_Rs))])
	      weights_Rs_sum = np.array([mean(a) for a in zip(*(weights_Rs))])
	      weighted_mean_Rs = [x for x in weighted_mean_Rs if x not in weighted_mean_Rs[i]] # Removing the contribution by the comparison star
	      weighted_mean_Rs_sum = np.array([np.median(a) for a in zip(*(weighted_mean_Rs))])
	      norm_ref_star_flux_weight_sum_mean_Rs.append(weighted_mean_Rs_sum / weights_Rs_sum)	  
	      REF_LCs.append(flux_r_norm[i]/norm_ref_star_flux_weight_sum_mean_Rs[i])
	      err_ref_Rs.append(np.sqrt(err_r_norm[i]**2 + np.median(sigma_weighted_mean_ref_stars)**2))	# IRAF errors
		  
	      REF_LCs[i] = REF_LCs[i]/np.median(REF_LCs[i])
 	  '''
 	  
 	  
	  ''' Finding a reference star of similar brightness for comparison sake '''
	  '''
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

	    for i in range(len(mag_r)):			# Finding the REF star closest in brightness.
		  if mag_refs[i] == min(mag_refs, key=lambda x:abs(x-np.median(flux_T_MAG))):
		    ref_number = i					# ref_number indicates the choosen ref star

	  print "\nTarget and Comparison star mag difference:\t",round(abs(np.median(flux_T_MAG) - mag_refs[ref_number]),2),'\t(Target Brighter)\n' if ((np.median(flux_T_MAG) - mag_refs[ref_number]) >= 0) else '\t(Ref Brighter)\n'
	  '''
