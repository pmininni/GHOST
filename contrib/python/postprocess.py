import numpy as np
import glob as glob
import string as string
import matplotlib.pyplot as plt

# Reads binary files in a directory and does
# postprocessing (e.g., to visualize later with
# VAPOR). Note that for big runs, GHOST can do some
# automatic postprocessing (e.g., compute vorticity)
# at run time.
# Execute in ipython with '%run postprocess.py'

# Path to the binary data
path = '../../3D/bin/outs/'

# Box size
Lx = 2*np.pi
Ly = 2*np.pi

# Spatial resolution
NX = 128
NY = 128
NZ = 128
dx = Lx/NX
dy = Ly/NX
shape = (NX,NY,NZ)

# Reads binary files, computes vertical vorticity
# using centered finite differences, and saves in 
# a new binary file named 'wz.NNNN.out'
filelist = sorted(glob.glob(path+'vx.*.out'))
for file in filelist:
  vx = np.fromfile(file,dtype=np.float32).reshape(shape,order='F')
  ind = string.lstrip(file,path+'vx.').rstrip('.out')
  str = 'vy.'+ind+'.out'
  vy = np.fromfile(path+str,dtype=np.float32).reshape(shape,order='F')
  adv = np.roll(vy,-1,axis=0)
  ret = np.roll(vy,1,axis=0)
  vy = (adv-ret)/(2*dx) # dv_y/dx
  adv = np.roll(vx,-1,axis=1)
  ret = np.roll(vx,1,axis=1)
  vx = (adv-ret)/(2*dy) # dv_x/dy 
  wz = vy-vx
  wz.tofile(path+'wz.'+ind+'.out')
