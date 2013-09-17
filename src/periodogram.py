import numpy as np
import scipy as sp
from scipy import signal
from scipy.signal import lombscargle
import matplotlib.pyplot as plt


def periodogram(x,y,start,stop,number):
  f = np.linspace(start,stop,number)
  normval = x.shape[0]
  pgram = sp.signal.lombscargle(x, y, f)
  max_value = np.sqrt(4*(pgram/normval)).max()
  return pgram,f,normval,max_value

def getSignificance(wk1, wk2, nout, ofac):
    """ returns the peak false alarm probabilities
    Hence the lower is the probability and the more significant is the peak
    """
    expy = exp(-wk2)                   
    effm = 2.0*(nout)/ofac              
    sig = effm*expy
    ind = (sig > 0.01).nonzero()
    sig[ind] = 1.0-(1.0-expy[ind])**effm
    return sig

#periodogram(np.array([0,2,3]),np.array([2,5,7]),0.,10.,100.0)



fs = 10e3
N = 1e5
amp = 2*np.sqrt(2)
freq = 1234.0
noise_power = 0.001 * fs / 2
time = np.arange(N) / fs
x = amp*np.sin(2*np.pi*freq*time)
x += np.random.normal(scale=np.sqrt(noise_power), size=time.shape)
print x
'''
pgram = sp.signal.lombscargle(x, y, f)

plt.subplot(2, 1, 1)
plt.plot(x, y, 'b+')


plt.subplot(2, 1, 2)
plt.plot(f, np.sqrt(4*(pgram/normval)))
plt.show()
'''