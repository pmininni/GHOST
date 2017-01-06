import numpy as np
import matplotlib.pyplot as plt

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
path = '../../3D/bin/outs/'

# Spatial resolution
NX = 128
NY = 128
NZ = 128
shape = (NX,NY,NZ)

# Reads binary files
vx = np.fromfile(path+'vx.0001.out',dtype=np.float32).reshape(shape,order='F')

# Show a horizontal cut of the field in the middle of the box
plt.figure(1)
plt.imshow(vx[:,:,N/2])
plt.xlabel('x')
plt.ylabel('y')
plt.show()
