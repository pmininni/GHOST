!=================================================================
! MODULES for POSIX binary I/O
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!=================================================================

  MODULE iovar
!
! Change the length of the string 'node' to change the number 
! of characters used to write the number of the processor.
! The format fmtnod should be consistent with the length of 
! the string, e.g. if len=5 then fmtnod = '(i5.5)'.
      TYPE IOPLAN
         INTEGER               :: n,jsta,jend
         CHARACTER(len=3)      :: node
      END TYPE IOPLAN
      CHARACTER(len=6)         :: fmtnod = '(i3.3)'
      SAVE

  END MODULE iovar
!=================================================================
