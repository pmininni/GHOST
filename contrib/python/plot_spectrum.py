import numpy as np
import glob as glob
import matplotlib.pyplot as plt

# Plots energy spectra as a function of time
# Execute in ipython with '%run plot_spectrum.py'

# Path to the data
path = '../../3D/bin/'

# Reads and plots all spectra in the directory.
# We only plot one every five spectra, starting
# from the second.
filelist = sorted(glob.glob(path+'kspectrum.*.txt'))
nfiles = np.size(filelist)
plt.figure(1)
for i in range(1, nfiles, 5):
  ene = np.loadtxt(filelist[i])
  plt.loglog(ene[:,0],ene[:,1])
nmax = 2*(np.size(ene[:,0])-1)
plt.xlim(1,nmax/3)
plt.xlabel('k')
plt.ylabel('E(k)')
plt.show()
