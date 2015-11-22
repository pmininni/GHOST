PRO lapla_findiff,in,out,dx,DY=dy,DZ=dz
;+
; NAME:
;       LAPLA_FINDIFF
;
; PURPOSE:
;       Computes the laplacian using sixth order finite 
;       differences in regular Cartesian grids
;
; CALLING SEQUENCE:
;       LAPLA_FINDIFF,IN,OUT,DX
;
; PARAMETERS:
;       IN[in]:   3D array with the field component
;       OUT[out]: 3D array with the laplacian
;       DX[in]:   spatial step of the grid in the direction DIR
;
; KEYWORD PARAMETERS:
;       DY[in]:    spatial step of the grid in the y direction.
;                  If not set, its value is set to the value of DX
;       DZ[in]:    spatial step of the grid in the z direction.
;                  If not set, its value is set to the value of DX
;
; COMMON BLOCKS:
;       None
;
; AUTHOR:
;       Pablo Daniel Mininni
;-

on_error,2                      ;return to caller if an error occurs
if keyword_set(dy) eq 0 then dy=dx
if keyword_set(dz) eq 0 then dz=dx

s = size(in)                    ;size of the input array
out = in                        ;output has the same size than input
d2x = out & d2y = out           ;temporary arrays

dir = 1                         ;second derivative in the x direction
for i = 0,2 do begin            ;forward differences near the boundary
    d2x(i,*,*) = (2*in(i,*,*)-5*in(i+1,*,*)+4*in(i+2,*,*)-in(i+3,*,*))/dx^2
endfor
for i = 3,s(dir)-4 do begin     ;centered differences
    d2x(i,*,*) = (2*in(i-3,*,*)-27*in(i-2,*,*)+270*in(i-1,*,*)-490*in(i,*,*) $
                 +270*in(i+1,*,*)-27*in(i+2,*,*)+2*in(i+3,*,*))/(180*dx^2)
endfor
for i = s(dir)-3,s(dir)-1 do begin ;backward differences near the boundary
    d2x(i,*,*) = (-in(i-3,*,*)+4*in(i-2,*,*)-5*in(i-1,*,*)+2*in(i,*,*))/dx^2
endfor

dir = 2                         ;second derivative in the y direction
for i = 0,2 do begin            ;forward differences near the boundary
    d2y(*,i,*) = (2*in(*,i,*)-5*in(*,i+1,*)+4*in(*,i+2,*)-in(*,i+3,*))/dy^2
endfor
for i = 3,s(dir)-4 do begin     ;centered differences
    d2y(*,i,*) = (2*in(*,i-3,*)-27*in(*,i-2,*)+270*in(*,i-1,*)-490*in(*,i,*) $
                 +270*in(*,i+1,*)-27*in(*,i+2,*)+2*in(*,i+3,*))/(180*dy^2)
endfor
for i = s(dir)-3,s(dir)-1 do begin ;backward differences near the boundary
    d2y(*,i,*) = (-in(*,i-3,*)+4*in(*,i-2,*)-5*in(*,i-1,*)+2*in(*,i,*))/dy^2
endfor

dir = 3                         ;second derivative in the z direction
for i = 0,2 do begin            ;forward differences near the boundary
    out(*,*,i) = (2*in(*,*,i)-5*in(*,*,i+1)+4*in(*,*,i+2)-in(*,*,i+3))/dz^2
endfor
for i = 3,s(dir)-4 do begin     ;centered differences
    out(*,*,i) = (2*in(*,*,i-3)-27*in(*,*,i-2)+270*in(*,*,i-1)-490*in(*,*,i) $
                 +270*in(*,*,i+1)-27*in(*,*,i+2)+2*in(*,*,i+3))/(180*dz^2)
endfor
for i = s(dir)-3,s(dir)-1 do begin ;backward differences near the boundary
    out(*,*,i) = (-in(*,*,i-3)+4*in(*,*,i-2)-5*in(*,*,i-1)+2*in(*,*,i))/dz^2
endfor

out = out+d2x+d2y

end
