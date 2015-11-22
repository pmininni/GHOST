; Plots energy spectra as a function of time
; Type '.r plot_spectrum' in IDL to execute

; Path to the data and number of files to read
path = '../../3D/bin/'
nfiles = 20

; Spatial resolution
N = 128

; Array with wavenumbers starting at 1
k = findgen(N/2+1)+1

; Reads and plots all spectra in the directory.
; We only plot one every five spectra, starting
; from the second.
FORMAT = '(%"kspectrum.%4.4d.txt")'
window,0
for i = 2,nfiles,5 do begin
  loadtxt,path+string(i,FORMAT=FORMAT),ene
  if (i eq 2) then begin
     plot,k,ene,/xlog,/ylog,xr=[1,N/3],xs=1,xtitle='k',ytitle='E(k)'
  endif else begin
     oplot,k,ene
  endelse
endfor

end
