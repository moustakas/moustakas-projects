;+
; NAME:
;   TRGB_LFMODEL  -- likelihood distribution function, smoothed
; PURPOSE:
;   function to evaluate distribution function consisting of 3
;   linear pieces. mc is the cutoff magnitude, input, relative to
;   minmag specified at edge of box. cutfactor in range 0-1 is the
;   range of the sharp drop at TRGB. beta [optional] is slope of
;   bright  end of LF.
;   This routine is called by TRGB_MAXLIKE
; CALLING SEQUENCE:
;   F = TRGB_LFMODEL(error, mc, cutfactor,[alpha=alpha, beta=beta])
; INPUTS:
;   error-- vector that gives error (in mags) versus magnitude (500 elements
;           long, with 100 steps per magnitude)
;   mc -- magnitude of TRGB
;   cutfactor -- strength of discontinuity
; MODIFICATION HISTORY:
; md 22jun00, 15sep00
;    22sep00 -- add beta a additional fit parameter, rather than
;               default beta=.8
; bjm 29sep00
;-
function TRGB_LFMODEL_bm, error,  mc, cutfactor, eminus, eplus, nsmooth, a1, s1,$
                          beta = beta, alpha = alpha

nn = 500
dm = 100
g = FLTARR(nn)
mag = FINDGEN(nn)/dm ; magnitudes distributed in range minmag, minmag+5

IF NOT KEYWORD_SET(alpha) THEN  alpha = 0.30  ;optional input parameter
IF NOT KEYWORD_SET(beta) THEN  beta=0.8  ;optional input parameter

x = mag - mc  ; get offset
g = 10.^(alpha*x)*(x ge 0.) + 10.^(-cutfactor)*10.^(beta*x)*(x lt 0)

gs = FLTARR(500)
FOR i = 0, nn-1 DO gs[i]=TOTAL(g[eminus[i]:eplus[i]])/nsmooth[i] ;boxcar smooth

gs[a1[0]:a1[s1-1]] = g[a1[0]:a1[s1-1]]

stop
return,gs
end
