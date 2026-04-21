# GHOST

GHOST (the Geophysical High-Order Suite for Turbulence) is an accurate and highly scalable pseudospectral code that solves a variety of PDEs often encountered in studies of turbulent flows. It is mainly developed by Pablo Mininni and Duane Rosenberg, with contributions from many users. The code uses a hybrid parallelization method combining MPI and OpenMP (GHOST 2.x also has support for GPUs using CUDA and Intel MKL, while more general GPU support is planned for the current version). This allows the code to run efficiently in laptops and small clusters, as well as to scale up to over 100,000 CPU cores in production runs in supercomputers. The hybrid parallelization method used in GHOST was recognized by two awards: an NCAR/CISL award, and the best paper award (technology track) at the TeraGrid 2010 conference.

The most recent release (version 3) is object oriented and provides many new functionalities, while dropping support for some old libraries, compilers, and solvers. We also maintain stable branches with GHOST ver. 1 and 2, that provide more solvers, and GPU support in version 2.

GHOST can solve PDEs in periodic domains to tackle many problems including:

* Compressible and incompressible hydrodynamic equations.
* Compressible and incompressible magnetohydrodynamic (MHD) equations.
* Hall-MHD equations and other kinetic plasma effects.
* Rotating flows.
* Rotating and stratified (dry and moist) flows using the Boussinesq approximation.
* Passive scalars.
* Superfluids and Bose-Einstein condensates using the Gross-Pitaevskii Equation (GPE).
* Ginzburg-Landau equations.
* Integration of many types of particles.

GHOST is mostly written in Fortran 2008, with some C bindings. To build GHOST only a C and a Fortran compiler are required.  More details on how to compile GHOST can be found in src/README. GHOST also includes examples of scripts in Python, Matlab, IDL, and Mathematica, to do post-analysis of the numerical simulations. Moreover, it can be used together with VAPOR (the Visualization and Analysis Platform for Ocean, Atmosphere, and Solar Researchers) to do still frame renderings and animations.  VAPOR provides an interactive environment for 3D visualization of GHOST outputs, allowing powerful analysis of turbulent flows in systems with 3D graphics cards (https://www.vapor.ucar.edu/).
