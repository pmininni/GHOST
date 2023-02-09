import numpy as np

# Reads particles files in binary format. All particles
# files have the same format: First the total number of
# particles, then the time, and then an array with shape
# (Npart,3) with the positions (or velocities, depending
# of the file name) of all particles at the given time.
# Execute in ipython with '%run load_particles.py'

# Path to the binary data
path = '../../3D/bin/outs/'

# Reads positions of any type of particles
pos  = np.fromfile(path+'/xlg.00000001.lag',dtype=np.float32)
N    = pos[0]
time = pos[1]
pos  = pos[2:].reshape((int(N),3))

# We now have all the positions in the array pos. We can
# separate x, y, and z in different arrays
x = pos[:,0]
y = pos[:,1]
z = pos[:,2]

# Reads velocities of inertial particles
vel  = np.fromfile(path+'/xlg.00000001.lag',dtype=np.float32)
N    = vel[0]
time = vel[1]
vel  = vel[2:].reshape((int(N),3))

# We can do the same with the velocity and separate
# components in different arrays
vx = vel[:,0]
vy = vel[:,1]
vz = vel[:,2]
