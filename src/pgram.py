#Code for the one planet system using the optimize sub-routine. Much much faster than V1.0

from numpy import sin,cos,zeros,arctan2,arange
from scipy.optimize import curve_fit,fmin,fmin_tnc,fmin_l_bfgs_b,fmin_cobyla
from pylab import fromfile
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy
from scipy import stats
from scipy.stats import kstest, ks_2samp
import math
import numpy
from numpy.fft import *

'''
HD21693.txt
HD21749.txt
HD4308.txt
'''


z=fromfile('data/variables/35b.txt', dtype=float, sep=' ')
z_ref=fromfile('data/variables/35rb.txt', dtype=float, sep=' ')
l=len(z)/3
l_ref = len(z_ref)/3
z.shape=n=(len(z)/3),3
z_ref.shape=n_ref=(len(z_ref)/3),3
verr=z[0:l,2]*1000
y=z[0:l,1]*1000
x=z[0:l,0]+2450000
verr_ref=z_ref[0:l,2]*1000
y_ref=z_ref[0:l_ref,1]*1000
x_ref=z_ref[0:l_ref,0]+2450000
Ms=1.053*1.99E30
start=min(x)

ym=y.mean()
i=0
while i<len(x):
  i+=1
  

#------------------------------------------------------------------------------  
#-----------------------------------FUNCTIONS----------------------------------
def __spread__(y, yy, n, x, m):
  """
  Given an array yy(0:n-1), extirpolate (spread) a value y into
  m actual array elements that best approximate the "fictional"
  (i.e., possible noninteger) array element number x. The weights
  used are coefficients of the Lagrange interpolating polynomial
  Arguments:
    y : 
    yy : 
    n : 
    x : 
    m : 
  Returns:
  """
  nfac=[0,1,1,2,6,24,120,720,5040,40320,362880]
  if m > 10. :
    print 'factorial table too small in spread'
    return
  ix=long(x)
  if x == float(ix): 
    yy[ix]=yy[ix]+y
  else:
    ilo = long(x-0.5*float(m)+1.0)
    ilo = min( max( ilo , 1 ), n-m+1 ) 
    ihi = ilo+m-1
    nden = nfac[m]
    fac=x-ilo
    for j in range(ilo+1,ihi+1): fac = fac*(x-j)
    yy[ihi] = yy[ihi] + y*fac/(nden*(x-ihi))
    for j in range(ihi-1,ilo-1,-1):
      nden=(nden/(j+1-ilo))*(j-ihi)
      yy[j] = yy[j] + y*fac/(nden*(x-j))

def fasper(x,y,ofac,hifac, MACC=4):
  """ function fasper
    Given abscissas x (which need not be equally spaced) and ordinates
    y, and given a desired oversampling factor ofac (a typical value
    being 4 or larger). this routine creates an array wk1 with a
    sequence of nout increasing frequencies (not angular frequencies)
    up to hifac times the "average" Nyquist frequency, and creates
    an array wk2 with the values of the Lomb normalized periodogram at
    those frequencies. The arrays x and y are not altered. This
    routine also returns jmax such that wk2(jmax) is the maximum
    element in wk2, and prob, an estimate of the significance of that
    maximum against the hypothesis of random noise. A small value of prob
    indicates that a significant periodic signal is present.  
  Reference: 
    Press, W. H. & Rybicki, G. B. 1989
    ApJ vol. 338, p. 277-280.
    Fast algorithm for spectral analysis of unevenly sampled data
    (1989ApJ...338..277P)
  Arguments:
      X   : Abscissas array, (e.g. an array of times).
      Y   : Ordinates array, (e.g. corresponding counts).
      Ofac : Oversampling factor.
      Hifac : Hifac * "average" Nyquist frequency = highest frequency
           for which values of the Lomb normalized periodogram will
           be calculated.      
   Returns:
      Wk1 : An array of Lomb periodogram frequencies.
      Wk2 : An array of corresponding values of the Lomb periodogram.
      Nout : Wk1 & Wk2 dimensions (number of calculated frequencies)
      Jmax : The array index corresponding to the MAX( Wk2 ).
      Prob : False Alarm Probability of the largest Periodogram value
      MACC : Number of interpolation points per 1/4 cycle
            of highest frequency
  History:
    02/23/2009, v1.0, MF
      Translation of IDL code (orig. Numerical recipies)
  """
  #Check dimensions of input arrays
  n = long(len(x))
  if n != len(y):
    print 'Incompatible arrays.'
    return
  nout  = 0.5*ofac*hifac*n
  nfreqt = long(ofac*hifac*n*MACC)   #Size the FFT as next power
  nfreq = 64L             # of 2 above nfreqt.
  while nfreq < nfreqt: 
    nfreq = 2*nfreq
  ndim = long(2*nfreq)  
  #Compute the mean, variance
  ave = y.mean()
  ##sample variance because the divisor is N-1
  var = ((y-y.mean())**2).sum()/(len(y)-1) 
  # and range of the data.
  xmin = x.min()
  xmax = x.max()
  xdif = xmax-xmin

  #extirpolate the data into the workspaces
  wk1 = zeros(ndim, dtype='complex')
  wk2 = zeros(ndim, dtype='complex')

  fac  = ndim/(xdif*ofac)
  fndim = ndim
  ck  = ((x-xmin)*fac) % fndim
  ckk  = (2.0*ck) % fndim

  for j in range(0L, n):
    __spread__(y[j]-ave,wk1,ndim,ck[j],MACC)
    __spread__(1.0,wk2,ndim,ckk[j],MACC)

  #Take the Fast Fourier Transforms
  wk1 = ifft( wk1 )*len(wk1)
  wk2 = ifft( wk2 )*len(wk1)

  wk1 = wk1[1:nout+1]
  wk2 = wk2[1:nout+1]
  rwk1 = wk1.real
  iwk1 = wk1.imag
  rwk2 = wk2.real
  iwk2 = wk2.imag
  
  df  = 1.0/(xdif*ofac)
  
  #Compute the Lomb value for each frequency
  hypo2 = 2.0 * abs( wk2 )
  hc2wt = rwk2/hypo2
  hs2wt = iwk2/hypo2

  cwt  = numpy.sqrt(0.5+hc2wt)
  swt  = numpy.sign(hs2wt)*(numpy.sqrt(0.5-hc2wt))
  den  = 0.5*n+hc2wt*rwk2+hs2wt*iwk2
  cterm = (cwt*rwk1+swt*iwk1)**2./den
  sterm = (cwt*iwk1-swt*rwk1)**2./(n-den)

  wk1 = df*(arange(nout, dtype='float')+1.)
  wk2 = (cterm+sterm)/(2.0*var)
  pmax = wk2.max()
  jmax = wk2.argmax()


  #Significance estimation
  #expy = exp(-wk2)          
  #effm = 2.0*(nout)/ofac       
  #sig = effm*expy
  #ind = (sig > 0.01).nonzero()
  #sig[ind] = 1.0-(1.0-expy[ind])**effm

  #Estimate significance of largest peak value
  expy = numpy.exp(-pmax)          
  effm = 2.0*(nout)/ofac       
  prob = effm*expy

  if prob > 0.01: 
    prob = 1.0-(1.0-expy)**effm

  return wk1,wk2,nout,jmax,prob

def getSignificance(wk1, wk2, nout, ofac):
  """ returns the peak false alarm probabilities
  Hence the lower is the probability and the more significant is the peak
  """
  expy = numpy.exp(-wk2)          
  effm = 2.0*(nout)/ofac       
  sig = effm*expy
  ind = (sig > 0.01).nonzero()
  sig[ind] = 1.0-(1.0-expy[ind])**effm
  return sig
#------------------------------------------------------------------------------

def Tau(t,w):
  k1=0.0
  k2=0.0
  Ssin=0.0
  Scos=0.0
  while k1<len(t):
    Ssin=Ssin+sin(2*w*t[k1])
    k1=k1+1
  while k2<len(t):
    Scos=Scos+cos(2*w*t[k2])
    k2=k2+1
  return (math.atan(Ssin/Scos))/(2*w)

def sum1(x,t,w):
  sum=0
  i1=0
  while i1<len(t):
    sum=sum+x[i1]*cos(w*t[i1])
    i1=i1+1
  sum=sum**2
  return sum
 
def sum2(t,w):
  sum=0
  i2=0
  while i2<len(t):
    sum=sum+(cos(w*t[i2]))**2
    i2=i2+1
  return sum

def sum3(x,t,w):
  sum=0
  i3=0
  while i3<len(t):
    sum=sum+x[i3]*sin(w*t[i3])
    i3=i3+1
  sum=sum**2
  return sum

def sum4(t,w):
  sum=0
  i4=0
  while i4<len(t):
    sum=sum+(sin(w*t[i4]))**2
    i4=i4+1
  return sum
    
def theta(x,e,P,To):
  #j=0
  #while j<len(x):
  #  x[j]=x[j]-start_time
  #  j+=1
  EA=zeros( (len(x)) )
  MA=zeros( (len(x)) )
  i=0
  while i<len(x):
    MA[i]=2*scipy.pi*(x[i]-To)/P
    EA[i]=eccan(e,MA[i])
    i+=1
  CosTheta=((scipy.cos(EA))-e)/(1.0-e*scipy.cos(EA))
  SinTheta=(scipy.sin(EA)*scipy.sqrt(1.0-e**2)/(1.0-e*scipy.cos(EA)))
  True_A=arctan2(SinTheta,CosTheta)
  return True_A

def pgram(x,y,max,step):
  #Create an array of values for the frequency. Using my method.
  #This value determines the range of frequencies (i.e. period range) to search.
  Ang_range=max
  Pmax=x.max()
  w=arange(0,Ang_range*math.pi/Pmax,step)
  Pw=zeros( (len(w),) )
  j=0
  #Build a LS periodogram using predefined summations,
  while j<len(w):
    tau=Tau(x,w[j])
    Pw[j]=0.5*(sum1(y,x-tau,w[j])/sum2(x-tau,w[j])+sum3(y,x-tau,w[j])/sum4(x-tau,w[j]))
    j=j+1
  #Pull out the maximum power so that can find the associated frequency.
  int=1
  peak=0.0
  freq=0.0
  while int<len(Pw):
    if Pw[int]>peak:
      peak=Pw[int]
      freq=w[int]
    int+=1
  P=2*math.pi/freq
  w=w/(2*math.pi)
  print 'The radial velocity curve vaires over a peiod of',P,'days.'
  #Plot the periodogram
  plt.figure(1)
  plt.plot(w,Pw,'b-', markersize=2)
  plt.ylabel('P(Frequency)')
  plt.xlabel('Frequency')
  return P


  
 
#------------------------------------------------------------------------------  
#-----------------------------------FUNCTIONS-END------------------------------
#------------------------------------------------------------------------------


print "Run periodogram, y, or assume a period previously calculated, n?"
while True:
  answer='y'#raw_input("Enter either 'yes' or 'no'.")
  if answer=='y':
    break
  if answer=='n':
    break 
if answer=='y':
  print "Use the Num Rec. code or mine (r or m)?"
  while True:
    answer2='r'#raw_input("Enter either 'r' or 'm'")
    if answer2=='r':
      break
    if answer2=='m':
      break 

if (answer=='y' and answer2=='r'):
  fx,fy, derp, jmax, prob = fasper(x,y, 40., 1.)
  fx_ref,fy_ref, derp_ref, jmax_ref, prob_ref = fasper(x_ref,y_ref, 40., 1.)
  figP=plt.figure(1)
  fp=1.0/fx
  fp_ref=1.0/fx_ref
  axP=figP.add_subplot(1,1,1)
  #axP.set_xscale('log')
  #plt.xlim(min(fp),3000)
  plt.xlim(0.02,0.5)
  plt.plot(fp,fy,'k-', markersize=2)
  plt.plot(fp_ref,fy_ref,'r--', markersize=2)
  plt.ylabel('Power')
  plt.xlabel('Period [Days]')
  P=1./fx[jmax]
  P_ref=1./fx_ref[jmax_ref]
  print 'The false alarm probability for the largest peak is:',prob
elif (answer=='y' and answer2=='m'):
  maxrange=0.5
  step=0.005
  #maxrange=2000000
  #step=0.001
  P=pgram(x,y,maxrange,step)
  P_ref=pgram(x_ref,y_ref,maxrange,step)
print 'The period is:',P,'days (be it calc or preset).'
print 'The period of the referene is:',P_ref,'days (be it calc or preset).'
plt.savefig('periodogram.pdf')
plt.show()
#raw_input('FIN')