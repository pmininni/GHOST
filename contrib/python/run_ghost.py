# This file provides an example on how to run GHOST from Python using
# the dynamic library libghost_python. This method is useful to 
# integrate GHOST with ML workflows.

import ctypes
import os

# --------------------------------------------------------------------
# 1. Load the Shared Library (.so extension in Linux, .dylib in MacOS)
# --------------------------------------------------------------------
lib_path = os.path.abspath("../../lib/libghost_python.so")

# CDLL loads the library. RTLD_GLOBAL is often needed for MPI symbols 
# to resolve correctly across different libraries.
ghost = ctypes.CDLL(lib_path, mode=ctypes.RTLD_GLOBAL)

# --------------------------------------------------------------------
# 2. Define Argument Types
# --------------------------------------------------------------------
# Fortran expects arguments by reference, but since we use VALUE 
# in the Fortran interface for steps, we can pass integers directly.
ghost.ghost_run.argtypes = [ctypes.c_int]
ghost.ghost_run.restype  = None

# --------------------------------------------------------------------
# 3. Run the Simulation
# --------------------------------------------------------------------
print("--- Python: Initializing GHOST ---")
ghost.ghost_init()

print("--- Python: Running 10 steps ---")
ghost.ghost_run(10)

print("--- Python: Running 10 more steps ---")
ghost.ghost_run(10)

print("--- Python: Finalizing GHOST ---")
ghost.ghost_finalize()

print("--- Python: Done ---")
