#! /usr/bin/env python

# Stat definitions

import numpy as np
import numpy.ma as ma

def chi2(obs,err):
  obs = obs/np.median(obs)
  core = ((obs - 1.) / err)**2.
  return core.sum()/len(obs-1)

def chi2_master_ref(obs,err,comp,comp_err):
  obs = obs/np.median(obs)
  core = ((obs - comp) / np.sqrt(err**2+comp_err**2))**2.
  return core.sum()/len(obs-1)

def robust(obs,err):
  core = abs((obs-np.median(obs))/err)
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

def coloumn_median(B,A):
  maxlen = max(len(x) for x in B)
  if len(A) == maxlen:
    ans = A
  else:
    C = np.array([l+[np.nan]*(maxlen-len(l)) for l in B])
    dat = np.ma.masked_array(C,np.isnan(C))
    ans = np.median(dat,axis=0)
  return ans

def coloumn_weighted_mean(B,err,N_ref_stars):
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
