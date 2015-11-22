pro imagesc, a, x, y, WINDOW_SCALE = window_scale, ASPECT = aspect, $
	INTERP = interp, TITLE = title, XTITLE = xtitle , YTITLE = ytitle, $
        COLOR = color, BACKGROUND = background, NOERASE = noerase, $
        XSTYLE = xstyle, YSTYLE = ystyle, COLORBAR = colorbar
;+
; NAME:
;	IMAGESC
;
; PURPOSE:
;	Generate an image with a comand similar to imagesc in MATLAB and 
;       to imshow in matplotlib.
;
; CALLING SEQUENCE:
;	IMAGESC, A
;
; INPUTS:
;	A:	The two-dimensional array to display.
;
;       X,Y:    Two vectors containing x and y coordinates (optional)
;
; KEYWORD PARAMETERS:
; WINDOW_SCALE:	Set this keyword to scale the window size to the image
;		size. Otherwise, the image size is scaled to the window
;		size. This keyword is ignored when the output device has
;		scalable pixels (e.g., PostScript).
;
;	ASPECT:	Set this keyword to retain the image's aspect ratio.
;		Square pixels are assumed.  If WINDOW_SCALE is set, the 
;		aspect ratio is automatically retained.
;
;	INTERP:	If this keyword is set, bilinear interpolation is used if 
;		the image is resized.
;
;       COLORBAR: Set this keyword to draw a colorbar.
;
;       BACKGROUND: Background color.
;
;       COLOR: Axis and titles color.
;
;       [XY]TITLE: x- and y-axis titles.
;
;       [XY]STYLE: Axis style; 4 means no axis.
;
;       NOERASE: Retain previous figure when drawing.
;
; OUTPUTS:
;	No explicit outputs.
;
; COMMON BLOCKS:
;	None.
;
; PROCEDURE:
;	If the device has scalable pixels, then the image is written over
;	the plot window.
;
; AUTHOR:
;       Pablo Daniel Mininni
;-

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz[0] lt 2 then message, 'Parameter not 2D'

	;set window used by imagesc

if keyword_set(x) eq 0 then x = findgen(sz[1])
if keyword_set(y) eq 0 then y = findgen(sz[2])
if keyword_set(color) eq 0 then color = 1
if keyword_set(background) eq 0 then background = !d.n_colors-1
if keyword_set(xstyle) eq 0 then xstyle = 1
if keyword_set(ystyle) eq 0 then ystyle = 1

contour,a,x,y,/nodata, xstyle=4, ystyle = 4
if keyword_set(noerase) eq 0 then begin
  erase, background             ;Set background color
endif

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
if keyword_set(colorbar) then begin
   border = (px[1]-px[0])/10
   px[1] = px[1]-1.75*border
endif
swx = px[1]-px[0]		;Size in X in device units
swy = py[1]-py[0]		;Size in Y
six = float(sz[1])		;Image sizes
siy = float(sz[2])
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

if (!d.flags and 1) ne 0 then begin	;Scalable pixels?
  if keyword_set(aspect) then begin	;Retain aspect ratio?
				;Adjust window size
	if f ge 1.0 then swy = swy / f else swx = swx * f
  endif
        scl = !d.x_px_cm/!d.x_ch_size   ;resample using pixels per cm
        tv,poly_2d(bytscl(a),$  ;Have to resample image
                [[0,0],[1./scl,0]], [[0,1./scl],[0,0]],$
                keyword_set(interp),scl*six,scl*siy), $
                px[0],py[0],xsize = swx, ysize = swy, /device

endif else begin	        ;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tvscl,a,px[0],py[0]	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
	tv,poly_2d(bytscl(a),$	;Have to resample image
		[[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
		keyword_set(interp),swx,swy), $
		px[0],py[0]
	endelse			;window_scale
  endelse			;scalable pixels

contour,a,x,y,/nodata,/noerase,xstyle=xstyle,ystyle=ystyle, $ ; plot the axis
           pos = [px[0],py[0], px[0]+swx-1,py[0]+swy],/dev, $
           COLOR = color, XTITLE = xtitle, YTITLE = ytitle , $
           TITLE = title ;Draw axis

if keyword_set(colorbar) then begin                           ; plot the bar
   swini = swx+border/5
   swend = swx+4*border/10
   xrmp = findgen(10)
   yrmp = (max(a)-min(a))*findgen(256)/255+min(a)
   ramp = fltarr(10,256)
   for i = 0,9 do ramp(i,*) = yrmp
   if (!d.flags and 1) ne 0 then begin     ;Scalable pixels?
      tvscl,ramp,px[0]+swini,py[0],xsize = swend-swini, ysize = swy, /device
   endif else begin                        ;Not scalable pixels    
      tv,poly_2d(bytscl(ramp), $           ;Have to resample image
                [[0,0],[10/(swend-swini),0]], [[0,256/swy],[0,0]],$
                keyword_set(interp),swend-swini,swy), $
                px[0]+swini,py[0],/device
   endelse
   contour,ramp,xrmp,yrmp,/nodata,/noerase,xs=4,ys=4, $
           pos = [px[0]+swini,py[0],px[0]+swend-1,py[0]+swy],/dev, $
           COLOR = color,xticks = 1,xtickname = [' ',' ']
   axis,/xaxis,xticks = 1,xtickname = [' ',' '],ticklen = 0,COLOR = color
   axis,yaxis = 0,yticks = 1,ytickname = [' ',' '],ticklen = 0,COLOR = color
   yrmp = (max(a)-min(a))*findgen(5)/4+min(a)
   axis,yaxis = 1,yticks = 4,ytickv = yrmp,ticklen = 0,COLOR = color
endif

return
end
