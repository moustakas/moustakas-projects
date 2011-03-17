;+
; NAME:   
;            edge_detect
; PURPOSE: 
;    To detect via maximum likelihood analysis the edge of the 
;    TRGB, based on 
; CATEGORY:
;    Analysis routine
; CALLING SEQUENCE:
;    edge_detect, maglist, lnlikelihood
; INPUTS: 
;    maglist -- list of magnitudes for input stars 
; OUTPUT:
;    lnlike -- ln of likelihood curve as a function of m_c, the cutoff
;                  magnitude
; KEYWORDS:
;    maxmag -- maximum magnitude for search
;    minmag -- minimum magnitude for search
; MODIFICATION HISTORY:
; md- 21Jun00
;-
pro EDGE_DETECT, objname, maglist, lnlike, error, $ 
       minmag=minmag, maxmag=maxmag, hst=hst, halo=halo, core=core, ccut=ccut

; first input data
if not Keyword_set(minmag) then minmag=20.
if not Keyword_set(maxmag) then maxmag=25.
 
bmag=fix((maglist-minmag)*100.) ;binned magnitudes in range minmag,maxmag
numberbins=100*(maxmag-minmag)
bmag = bmag[WHERE((bmag ge 0) and (bmag lt numberbins))] ; subset of stars in selected magnitude range
 
Ntotal=n_elements(bmag)
lnlike =fltarr(50,20,6)
print, Ntotal,' stars within range', minmag, ' to ',maxmag

mag = FINDGEN(500)*0.01 + minmag
error = ERR_FUNC(objname, mag, hst=hst, halo=halo, core=core)

for l=0,5 do begin   ;loop over beta
  beta=l/10. +.5     ;slope of the bright end of LF
  for k= 0, 19 do begin
     cutfactor = k/20. +.1
     for i = 0, 49 do begin       ;loop over cutoff magnitude
        mc = i/10.                ; magnitude cutoff in range minmag, minmag+5
        gg = form_ggs(error, mc, cutfactor,beta) 
               ;note that this produces 500 points
        aloggg = alog(gg)
      lnlike[i,k,l] = -Ntotal*alog(total(gg[0:numberbins-1])) + $
                     total(aloggg(bmag))
     endfor     
  endfor
endfor


lnlike=lnlike-max(lnlike)
return
end
