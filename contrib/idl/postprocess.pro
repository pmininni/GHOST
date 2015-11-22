; Reads binary files in a directory and does
; postprocessing (e.g., to visualize later with
; VAPOR). Note that for big runs, GHOST can do some
; automatic postprocessing (e.g., compute vorticity)
; at run time.
; Type '.r postprocess' in IDL to execute

; Path to the binary data and number of files to read
path = '../../3D/bin/outs/'
nfiles = 4

; Spatial resolution
N = 128
dx = 2*!PI/N

; Reads binary files, computes vertical vorticity
; using sixth-order finite differences, and saves in
; a new binary file named 'wz.NNNN.out'
FORMAT = '(%".%4.4d.out")'
vx = fltarr(N,N,N)
vy = fltarr(N,N,N)
wz = fltarr(N,N,N)
for i = 1,nfiles do begin
  openr,1,path+'vx'+string(i,FORMAT=FORMAT)
  readu,1,vx
  close,1  
  openr,1,path+'vy'+string(i,FORMAT=FORMAT)
  readu,1,vy
  close,1  
  deriv_findiff,vy,wz,1,dx
  deriv_findiff,vx,vy,2,dx
  wz = wz-vy
  openw,1,path+'wz'+string(i,FORMAT=FORMAT)
  writeu,1,wz
  close,1
endfor

end
