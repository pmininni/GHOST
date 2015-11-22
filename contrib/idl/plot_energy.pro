; Plots energy as a function of time
; Assumes balance.txt is the output of an HD run
; Type '.r plot_energy' in IDL to execute

; Path to the data
path = '../../3D/bin/'

; Reads balance.txt
;  balance(:,1) = time
;  balance(:,2) = energy (v^2)
;  balance(:,3) = enstrophy (w^2)
;  balance(:,4) = energy injection rate
loadtxt,path+'balance.txt',balance,4

; Plots energy vs. time in a window
window,0
plot,balance(0,*),balance(1,*),xtitle='time',ytitle='Energy'

; Plots to an EPS file
set_plot,'ps'
device,xs=12,ys=10,filename='figure.eps',bits_per_pixel=8,/encapsulated 
plot,balance(0,*),balance(1,*),xtitle='time',ytitle='Energy'
device,/close
set_plot,'x'

end
