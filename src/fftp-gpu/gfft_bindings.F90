!=================================================================
! MODULE gfft
!
! Fortran interfaces to the device FFT library (cuFFT or hipFFT)
! through ISO_C_BINDING, so that the same source serves both
! vendors and any Fortran compiler. The plans are batched
! many-transform plans created with the C-order dimensions; the
! data layout and the (absent) normalization are those of FFTW, so
! the transforms are interchangeable with the FFTW ones of fftp-3.
! Execution routines receive device addresses, obtained by the
! caller inside a "target data use_device_addr" region.
!
! 2026 Pablo D. Mininni
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================
#include "gfft_defs.h"

      MODULE gfft
      USE ISO_C_BINDING
      IMPLICIT NONE

      INTERFACE
        INTEGER(C_INT) FUNCTION gfft_plan_many(plan,rank,n,inembed,   &
                    istride,idist,onembed,ostride,odist,ftype,batch)  &
                    BIND(C,NAME=C_GFFT_PLAN_MANY)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T :: plan
          INTEGER(C_INT), VALUE :: rank,istride,idist,ostride,odist
          INTEGER(C_INT), VALUE :: ftype,batch
          INTEGER(C_INT), INTENT(IN) :: n(*),inembed(*),onembed(*)
        END FUNCTION gfft_plan_many
        INTEGER(C_INT) FUNCTION gfft_exec_r2c(plan,in,out)            &
                    BIND(C,NAME=C_GFFT_EXEC_R2C)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T, VALUE :: plan
          TYPE(C_PTR), VALUE :: in,out
        END FUNCTION gfft_exec_r2c
        INTEGER(C_INT) FUNCTION gfft_exec_c2r(plan,in,out)            &
                    BIND(C,NAME=C_GFFT_EXEC_C2R)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T, VALUE :: plan
          TYPE(C_PTR), VALUE :: in,out
        END FUNCTION gfft_exec_c2r
        INTEGER(C_INT) FUNCTION gfft_exec_c2c(plan,in,out,dir)        &
                    BIND(C,NAME=C_GFFT_EXEC_C2C)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T, VALUE :: plan
          TYPE(C_PTR), VALUE :: in,out
          INTEGER(C_INT), VALUE :: dir
        END FUNCTION gfft_exec_c2c
        INTEGER(C_INT) FUNCTION gfft_destroy(plan) BIND(C,NAME=C_GFFT_DESTROY)
          IMPORT :: C_PTR,C_INT
          GFFT_HANDLE_T, VALUE :: plan
        END FUNCTION gfft_destroy
      END INTERFACE

      END MODULE gfft
