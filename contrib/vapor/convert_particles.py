from netCDF4 import Dataset
import numpy as np
import glob as glob

# This example shows how to convert GHOST particle data into
# VAPOR netCDF DCP files (Data Collection Particles). It
# assumes GHOST particle files are in binary format.

# Path to the binary data and simulation information
path = 'data'
NX = 1024    # Box resolution
NY = 1024
NZ = 512
LX = 2       # Box size
LY = 2
LZ = 1
N  = 1000000 # Number of particles
dim = 3      # Spatial dimensions
iun = 0      # 0 if data is in grid units (ilgwrtunit = 0), 1 otherwise 

# Looks for files
files = sorted(glob.glob(path+'/xlg.*.lag'))
ind   = 0

# Loads particles from the next available file
def StepFile() -> None:
    global pos, vel, ind
    pos = np.fromfile(files[ind],dtype=np.float32)
    pos = pos[2:].reshape((N,dim))
    if iun < 1:
       pos[:,0] = pos[:,0]*LX/NX
       pos[:,1] = pos[:,1]*LY/NY
       if dim > 2: pos[:,2] = pos[:,2]*LZ/NZ
    nmb = files[ind].lstrip(path+'/xlg.').rstrip('.lag')
    vel = np.fromfile(path+'/vlg.'+nmb+'.lag',dtype=np.float32)
    vel = vel[2:].reshape((N,dim))
    # We add two particles to set the domain size in VAPOR
    pos = np.append(pos, [[0 ,0 ,0 ],[LX,LY,LZ]], axis = 0)
    vel = np.append(vel, [[0 ,0 ,0 ],[0 ,0 ,0 ]], axis = 0)
    ind = ind+1

# Wraps an array in an additional dimension
def AddTimeDim(data: np.ndarray) -> np.ndarray:
    return np.expand_dims(data, axis=0)

def WriteTimestep(ts: float, positionData: np.ndarray, **kwargs: np.ndarray) -> None:
    dataset = Dataset(f"particles_{ts:08}.nc", "w", format="NETCDF4")
    particleCount:int = positionData.shape[0]
    dataset.createDimension("P", particleCount) # The P dimension represents the number of particles at this timestep
    dataset.createDimension("T", None) # Time dimension
    dataset.createDimension("axis", 3) # Utility dimension for packing 3 components for 3D particles into a single variable
    # Time coordinate
    T = dataset.createVariable("T", "f8", ("T",))
    T.units = "seconds"
    T[:] = np.array([ts])
    # 3D vars can be packed in a single variable by adding the axis dimension
    Position = dataset.createVariable("Position", "f4", ("T", "P", "axis"), zlib=True)
    # positionData is 2D (numParticles * axis) whereas Position is 3D (time * numParticles * axis)
    Position[:] = AddTimeDim(positionData)
    # Save all remaining particle properties passed in to nc file
    for name, data in kwargs.items():
        var = dataset.createVariable(name, "f4", ("T", "P", "axis")[0:data.ndim + 1], zlib=True)
        var[:] = AddTimeDim(data)
    dataset.close()

for ts in range(np.size(files)):
    StepFile()
    # Compute magnitude of velocity for each particle
    speed = np.sqrt(np.sum(vel**2, axis=-1))
    # Since 3-component properties such as velocity are common for 3D particles,
    # DCP allows packing them as a 2D array of size N_particles by 3
    # vel is an array of size Nx3 and speed is an array of size N
    WriteTimestep(ts+1, pos, vel=vel, speed=speed)
