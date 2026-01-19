! User-defined forcing scheme.
! A forcing scheme can be implemented in this file, which
! is used by the code when 'rand=3' is set in 'parameter.dat'.
! This scheme is executed every 'fstep' time steps. See the 
! folder 'examples' for an example. If not needed, this file
! can be left empty.
! Forcing arrays are complex (in Fourier space) of size
! (n,n,ista:iend) and are called:
!       (fx,fy,fz)   for the velocity 
!       (mx,my,mz)   for the e.m.f. (magnetic field)
!       (fs,fs1,...) for scalar fields
!       (fre,fim)    for quantum solvers
! You can use temporary real arrays R1-R3 of size
! (1:n,1:n,ksta:kend) and temporary complex arrays C1-C8 of
! size (n,n,ista:iend) to do intermediate computations,
! and two real arrays Faux1 and Faux2 of size (10) to store
! information of the history of the forcing if needed.
