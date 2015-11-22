import numpy as np
import matplotlib.pyplot as plt

# Plots energy as a function of time
# Assumes balance.txt is the output of an HD run
# Execute in ipython with '%run plot_energy.py'

# Path to the data
path = '../../3D/bin/'

# Reads balance.txt
#  t   = time
#  ene = energy (v^2)
#  ens = enstrophy (w^2)
#  eps = energy injection rate
t, ene, ens, eps = np.loadtxt(path+'balance.txt',unpack=True)

# Plots energy vs. time in a window
plt.figure(1)
plt.plot(t,ene)
plt.xlabel('time')
plt.ylabel('Energy')
plt.show()

# Saves plot to an EPS file
plt.savefig('figure.eps', format='eps', dpi=1000)
