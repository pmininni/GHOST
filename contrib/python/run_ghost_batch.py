# This file provides an example on how to run GHOST from Python using
# the dynamic library libghost_python. This method is useful, e.g., to 
# integrate GHOST within ML workflows.
#
# This script assumes a GHOST input file ('parameter.inp') is placed 
# in the same directory where python is called. To run this script in
# parallel as a non-interactive job, e.g., with 4 MPI tasks, do
# 'mpirun -np 4 python3 run_ghost.py'

import numpy as np
import ctypes
import os
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize   = False
from mpi4py import MPI
import platform
system = platform.system()
if   system == "Darwin":
    ext = ".dylib"
elif system == "Windows":
    ext = ".dll"
else:
    ext = ".so"

# --------------------------------------------------------------------
# 1. Load the Shared Library (.so extension in Linux, .dylib in MacOS)
# --------------------------------------------------------------------
lib_path = os.path.abspath(f"../../lib/libghost_python{ext}")

# CDLL loads the library. RTLD_GLOBAL is often needed for MPI symbols 
# to resolve correctly across different libraries.
ghost = ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)

# --------------------------------------------------------------------
# 2. Define Argument Types
# --------------------------------------------------------------------
# Fortran expects arguments by reference, but since we use VALUE 
# in the Fortran interfaces, we can pass arguments directly.
ghost.ghost_init.argtypes           = [ctypes.c_char_p]
ghost.ghost_init.restype            = None
ghost.ghost_run.argtypes            = [ctypes.c_int]
ghost.ghost_run.restype             = None
ghost.ghost_set_dt.argtypes         = [ctypes.c_double]
ghost.ghost_set_dt.restype          = None
ghost.ghost_get_real_size.argtypes  = [ctypes.c_int]
ghost.ghost_get_real_size.restype   = ctypes.c_int
ghost.ghost_get_real_field.argtypes = [ctypes.c_int]
ghost.ghost_get_real_field.restype  = ctypes.c_void_p

# --------------------------------------------------------------------
# 3. Run the simulation
# --------------------------------------------------------------------
print("--- Python: Initializing GHOST ---")
ghost.ghost_init('parameter.inp'.encode()) # We pass the string by bytes
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print("--- Python: Running 10 steps ---")
ghost.ghost_run(10)

# --------------------------------------------------------------------
# 4. Do some data analysis
# --------------------------------------------------------------------
nx  = ghost.ghost_get_real_size(1)
ny  = ghost.ghost_get_real_size(2)
nz  = ghost.ghost_get_real_size(3)
ptr = ghost.ghost_get_real_field(1)
if rank == 0:
    import matplotlib.pyplot as plt
    # Convert raw pointer to NumPy view. We assume GHOST was compiled with
    # single precision. When using double precision, we must use c_double.
    array_type = ctypes.c_float * (nx * ny * nz)
    buf        = array_type.from_address(ptr)
    field      = np.ctypeslib.as_array(buf).reshape((nx, ny, nz), order='F')
    plt.imshow(field[:,:,0])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig("figure.pdf")

# --------------------------------------------------------------------
# 4. Continue the simulation with a different time step
# --------------------------------------------------------------------
print("--- Python: Running 10 more steps ---")
ghost.ghost_set_dt(1e-2)
ghost.ghost_run(10)

print("--- Python: Finalizing GHOST ---")
ghost.ghost_finalize()

print("--- Python: Done ---")
