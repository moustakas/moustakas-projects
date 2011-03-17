function form_gg, mc, cutfactor
; silly function to evaluate distribution function consisting of 3
; linear pieces
; mc is the cutoff magnitude, input, relative to minmag specified at
; edge of box
; cutfactor in range 0-1 is the range of the sharp drop at TRGB.
; called by edge_detect and by findit
;
; md 22jun00
;
;bmc=(mc-20.)*10.
;if(bmc gt 49 or bmc lt 0) then  print, mc, ' magnitude cut out of range'
g= fltarr( 500)
mag=findgen(500)/100.  ; magnitudes distributed in range minmag, minmag+5
alpha=0.33  ; fiddle with these later
;beta=.8
beta=.8
gamma= 5.
x=mag-mc  ; get offset
g= 10.^(alpha*x)*(x ge 0.) + 10.^(gamma*x)*(x lt 0 and x ge -cutfactor/gamma) $
   + 10.^(-cutfactor)*10.^(beta*(x+cutfactor/gamma))*(x lt -cutfactor/gamma)

return,g
end



