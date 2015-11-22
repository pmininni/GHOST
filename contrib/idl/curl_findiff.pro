PRO curl_findiff,inx,iny,inz,outx,outy,outz,dx,DY=dy,DZ=dz
;+
; NAME:
;       CURL_FINDIFF
;
; PURPOSE:
;       Computes the curl of a vector field using sixth 
;       order finite differences in regular Cartesian grids
;
; CALLING SEQUENCE:
;       CURL_FINDIFF,INX,INY,INZ,OUTX,OUTY,OUTZ,DX
;
; PARAMETERS:
;       INX[in]:   3D array with the x component of the field
;       INY[in]:   3D array with the y component of the field
;       INZ[in]:   3D array with the z component of the field
;       OUTX[out]: 3D array with the x component of the curl
;       OUTY[out]: 3D array with the y component of the curl
;       OUTZ[out]: 3D array with the z component of the curl
;       DX[in]:    spatial step of the grid in the x direction
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

deriv_findiff,inz,aux1,2,dz       ;x component of the curl
deriv_findiff,iny,aux2,3,dy
outx = aux1-aux2

deriv_findiff,inx,aux1,3,dx       ;y component of the curl
deriv_findiff,inz,aux2,1,dz
outy = aux1-aux2

deriv_findiff,iny,aux1,1,dy       ;z component of the curl
deriv_findiff,inx,aux2,2,dx
outz = aux1-aux2

end
