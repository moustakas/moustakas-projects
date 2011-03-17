;+
; NAME:
;   form_ggs  -- likelihood distribution function, smoothed
; PURPOSE:
; silly function to evaluate distribution function consisting of 3
; linear pieces
; mc is the cutoff magnitude, input, relative to minmag specified at
; edge of box
; cutfactor in range 0-1 is the range of the sharp drop at TRGB.
; beta [optional] is slope of bright end of LF.
; This routine is called by edge_detect and by findit
; CALLING SEQUENCE:
;   f=form_ggs(error,mc,cutfactor,[beta])
; INPUTS:
;   error-- vector that gives error (in mags) versus magnitude (500 elements
;           long, with 100 steps per magnitude)
;   mc -- magnitude of TRGB
;   cutfactor -- strength of discontinuity
; MODIFICATION HISTORY:
;;
; md 22jun00, 15sep00
;    22sep00 -- add beta a additional fit parameter, rather than
;               default beta=.8
;-
function form_ggs, error,  mc, cutfactor,beta

nn=500
dm=100
g= fltarr(nn)
mag=findgen(nn)/dm ; magnitudes distributed in range minmag, minmag+5
alpha=0.30  ; fiddle with these later
;beta=.8
if (n_params() lt 4) then  beta=0.8  ;optional input parameter

x=mag-mc  ; get offset
g= 10.^(alpha*x)*(x ge 0.) + 10.^(-cutfactor)*10.^(beta*x)*(x lt 0)

ewidth=fix(error*dm*1.44)  ;cube root of 3 to get correct width sigma
eplus=(indgen(nn) +ewidth) < (nn-1)
eminus=(indgen(nn) -ewidth) > 0
nsmooth=eplus-eminus +1

gs=fltarr(500)
for i=0, nn-1 do gs[i]=total(g[eminus[i]:eplus[i]])/nsmooth[i] ;boxcar smooth

return,gs
end
