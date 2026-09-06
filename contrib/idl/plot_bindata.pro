; Reads a binary file and plots a cut in the x-y plane.
; Assumes we are reading a cubic box with resolution N^3.
; Type '.r plot_bindata' in IDL to execute

; Path to the binary data
path = '../../bin/outs/'

; Spatial resolution
N = 128

; Reads binary files
vx = fltarr(N,N,N)
openr,1,path+'vx.0001.out'
readu,1,vx
close,1

; Show a horizontal cut of the field in the middle of the box
window,0
imagesc,vx(*,*,N/2),xtitle='x',ytitle='y'

end
